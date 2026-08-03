# Conversão Rcpp — Modelo Poisson Polinomial de 2ª Ordem (DLM)

Este documento resume a conversão de R puro para C++/Rcpp dos cinco
amostradores MCMC do projeto, os padrões adotados, e o que falta para a
futura modularização em pacote R.

## Modelo

```
y_t        ~ Poisson(exp{theta_t1})
theta_t1   = theta_{t-1,1} + theta_{t-1,2} + omega_t1,   omega_t1 ~ N(0, W1)
theta_t2   =                 theta_{t-1,1} + omega_t2,   omega_t2 ~ N(0, W2)
theta_01   ~ N(mu_01, sigma2_01)
theta_02   ~ N(mu_02, sigma2_02)
1/W1 = phi1 ~ Gamma(nu_01, eta_01)
1/W2 = phi2 ~ Gamma(nu_02, eta_02)
```

## Arquivos

### Compartilhados (obrigatórios para todos os `.cpp`)

| Arquivo | Conteúdo |
|---|---|
| `utils.h` | Funções e estruturas comuns a **todos** os amostradores (não só Particle Gibbs) — ver referência completa abaixo. |

Todo `.cpp` faz `#include "utils.h"`; os dois arquivos precisam estar no
mesmo diretório (`sourceCpp()` resolve o include local automaticamente).

### Amostradores convertidos

| Algoritmo | `.cpp` | Driver R | Estratégia para $\theta_{t1}$ | Estratégia para $(\theta_{02},\theta_2)$ | Estratégia para $W_1, W_2$ |
|---|---|---|---|---|---|
| **PG-AS** | `pg_as.cpp` | `run_pg_as.R` | Bootstrap PF + Ancestral Sampling, *backward tracking* via genealogia | Chan (bloco estendido, exato) | Conjugado (Gibbs) |
| **AMH-Montoril** | `amh_montoril.cpp` | `run_amh_montoril.R` | Metropolis-within-Gibbs componente-a-componente, adaptativo (Robbins-Monro) | Chan (bloco estendido, exato) | Conjugado (Gibbs) |
| **SIR-Laplace** | `sir_laplace.cpp` | `run_sir_laplace.R` | IRLS/Laplace + Importance Sampling (ou passagem direta se $M_{is}=1$) | Chan (bloco estendido, exato) | Conjugado (Gibbs) |
| **SIR-Collapsed** | `sir_collapsed.cpp` | `run_sir_collapsed.R` | IRLS/Laplace + SIR | Chan (bloco estendido, exato) | **Colapsado**: MH independente com proposta Gamma calibrada por entropia cruzada; $\phi_1$ via verossimilhança IS, $\phi_2$ via verossimilhança exata |

`pg_apfB.cpp` / `run_pg_apfB.R` (dos primeiros turnos desta conversa) estão
**obsoletos** — o projeto migrou de PG-APF para PG-AS. Mantive nos anexos
apenas por histórico; não fazem parte do conjunto final.

## `utils.h` — referência de funções

```
logsumexp(x)                                    log-soma-exp de um vetor
log_p_yt(yt, theta)                             log-verossimilhança Poisson (com guarda de -Inf)
log_dnorm(x, mean, sd)                          log-densidade Normal

gibbs_sample_theta01(mu_01, sigma2_01,
                      theta_11, theta_02, W1)   passo conjugado de theta_01
gibbs_sample_phi1(nu_01, eta_01, theta_01,
                   theta1, theta_02, theta2, Tt) passo conjugado de phi1 (retorna PRECISÃO)
gibbs_sample_phi2(nu_02, eta_02, theta_02,
                   theta2, Tt)                  passo conjugado de phi2 (retorna PRECISÃO)

logpost_theta_t1(...) / logpost_theta_T1(...)   log-posteriori completa de theta_t1 (Montoril)
sample_theta_t1_mh(...) -> struct MHStep         passo de Metropolis (random-walk) para theta_t1

sample_index_from_logw(logw, K, w_buf)          1 índice via inversão de CDF (ordem natural)
sample_indices_from_logw(logw, K, w_buf,
                          cum_buf, out)         reamostragem SISTEMÁTICA (Kitagawa 1996) p/ K-1 partículas

struct TriLDLT { L, D }
tri_ldlt_factor(d, e)                           fatoração LDL^T tridiagonal (método de Chan)
tri_solve_A(fac, b, x)                          resolve A·x = b
tri_solve_Lt(fac, w, x)                         resolve L^T·x = w (usado para amostrar da Normal)
```

Todas as funções espelham nomes de `utils.R` 1:1 (quando existem no lado
R) — facilita o cross-reference na futura modularização.

## Padrões estabelecidos (aplicar em qualquer conversão futura)

**C++ / Rcpp**
- Parâmetros de estado inicial nomeados diretamente (`theta1`, `theta2`,
  `theta_01`, `theta_02`, `W1`, `W2`, ...), sem sufixo `_init` — são cópias
  locais por valor, mutadas em ordem.
- Vetores (`NumericVector` → `std::vector<double>`) exigem um bloco `{ }`
  aninhado para reusar o mesmo nome (parâmetro e variável local não podem
  coexistir no mesmo escopo — não é *shadowing*, é redeclaração ilegal em
  C++; um bloco aninhado cria escopo novo onde isso é permitido).
- `NumericMatrix`/`NumericVector` de histórico ficam declarados **fora**
  desse bloco (precisam sobreviver até o `return` final).
- Reamostragem: **sistemática** (Kitagawa 1996) quando são extraídos $K>1$
  índices de uma vez; **CDF-inversão de ordem natural** (`sample_index_
  from_logw`) para um único índice — nunca o `sample()`/`sample.int()` do
  R, que reordena as probabilidades internamente antes de acumular (ver
  nota de RNG abaixo).
- `RNGScope` sincroniza com `set.seed()` do R automaticamente — nenhuma
  semente precisa ser passada como parâmetro.
- `R::rgamma(shape, scale)` e `R::dgamma(x, shape, scale, log)` usam
  **escala**, não taxa — sempre passar `1.0/rate`.

**R (scripts `run_*.R`)**
- `setwd(dirname(this.path::this.path()))`.
- Bloco de carregamento de dados padronizado (ver qualquer `run_*.R` para
  o template exato — grade de cenários, seed derivada, `t_obs`).
- Um argumento por linha em chamadas de função; sem `;` para múltiplos
  comandos; sem alinhamento artificial de `<-`.
- Todos os gráficos/diagnósticos do script original preservados.

**Nota de RNG (importante, já mordeu duas vezes)**
`sample()`/`sample.int(n, size, prob=w)` do R **ordena `w` de forma
decrescente antes de acumular a CDF** (`ProbSampleReplace`, código-fonte
do R). Para bater bit-a-bit com uma implementação C++ que usa varredura em
ordem natural, é preciso *também* corrigir o lado R (adicionar uma função
`sample_index_from_logw()` equivalente) — feito em `sir_laplace.R` e
`sir_collapsed.R`. Sem essa correção, as duas implementações convergem
para a mesma posterior mas **não** batem trajetória-a-trajetória (validado
estatisticamente nesses casos antes da correção).

## Validação

Todos os cinco amostradores foram validados numericamente contra uma
referência fiel em R puro (ou contra o R original, quando disponível),
usando **a mesma seed**:

| Algoritmo | Resultado |
|---|---|
| PG-AS | ~$10^{-14}$ (bate com o `pg_as.R` autoritativo, sem nenhuma adaptação — já usava os mesmos algoritmos de amostragem) |
| AMH-Montoril | ~$10^{-16}$ a $10^{-14}$ |
| SIR-Laplace | ~$10^{-14}$ (após corrigir `sample()` → `sample_index_from_logw` no R) |
| SIR-Collapsed | ~$10^{-14}$, incluindo decisões de aceite/rejeite do MH **idênticas** e parâmetros da calibração CE idênticos |

Diferenças residuais na casa de $10^{-14}$–$10^{-18}$ são ruído de ponto
flutuante (ordem de operações distinta entre Cholesky esparsa do `Matrix`
e o solve tridiagonal direto), não divergência algorítmica.

## Desempenho (R puro vs. Rcpp, mesma seed, mesmos dados)

| Algoritmo | Speedup medido | Por quê |
|---|---|---|
| PG-AS | ~23–25× | Laço aninhado por partícula ($K\times T$ por iteração) — o que o R interpreta pior |
| SIR-Collapsed | ~20× | Maior volume de loops R aninhados por iteração (IRLS × MH-com-IS × SIR × verossimilhança exata) |
| AMH-Montoril | ~6,5–7,7× | Sem partículas; boa parte do tempo já ia para `Matrix::Cholesky` compilado |
| SIR-Laplace | ~5,2–7,4× | Idem — IRLS e $M_{is}$ pequenos, pouco código R "lento" a acelerar |

Padrão geral: quanto mais o algoritmo original dependia de laços R
aninhados por partícula/trajetória, maior o ganho; quanto mais o tempo já
ia para álgebra linear esparsa compilada, menor a margem — mas sempre uma
melhora líquida real, nunca marginal.

## Pendências para a modularização em pacote

Do `utils.R` original, ainda **não** foram portados para `utils.h` (porque
nenhum `.cpp` convertido até agora precisou deles diretamente — os
amostradores usam a máquina tridiagonal `TriLDLT` em vez da API baseada em
listas do R):

- `chan_build_chain`, `chan_build_static_objects[_ext]`,
  `make_chan_theta1_smoother_ext`, `make_chan_theta2_smoother_ext`,
  `chan_sample_from_build`, `chan_log_det_K0` — builders do método de Chan
  em R baseados em `Matrix::bandSparse`/`Cholesky`. Equivalente funcional
  já existe em C++ (`TriLDLT` + as montagens de $d,e,b$ específicas de
  cada bloco, inline em cada `.cpp`), mas a API não é 1:1 com o lado R.
- `cwmh_sample_theta1` (versão não-adaptativa do MH componente-a-componente
  — `sample_theta_t1_mh`/`logpost_theta_t1`/`logpost_theta_T1` já estão em
  `utils.h`, só falta a função de varredura completa sem a adaptação
  Robbins-Monro, caso algum outro algoritmo precise dela isoladamente).

Quando a modularização do pacote R começar, faz sentido decidir então se
vale a pena espelhar essa API baseada em listas em C++ (para paridade
literal com `utils.R`) ou se a abordagem `TriLDLT` genérica (mais enxuta,
já validada em 4 algoritmos diferentes) deve virar o padrão único do
pacote — na minha visão, a segunda opção é preferível, já que o `TriLDLT`
é estritamente mais simples e já provou ser reutilizável em blocos com
estruturas de diagonal bem diferentes (estendido vs. ancorado,
homocedástico vs. heterocedástico).
