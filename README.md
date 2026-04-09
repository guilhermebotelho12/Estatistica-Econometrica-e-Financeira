# Estatística Econométrica e Financeira - Trabalho Parte II (2024/25)

Este repositório contém a resolução e a implementação computacional da Parte II da unidade curricular de **Estatística Econométrica e Financeira** (Mestrado). [cite_start]O foco do trabalho é a estimação de parâmetros e a simulação de processos estocásticos baseados em Equações Diferenciais Estocásticas (EDE).

## Descrição do Problema
O estudo centra-se no processo $(X_t)_{t \ge 0}$ definido pela solução da EDE de Vasicek (Ornstein-Uhlenbeck):
$$dX_t = \alpha(\beta - X_t)dt + \sigma dW_t, \quad X_0 = x_0$$
Onde $\alpha$, $\beta$ e $\sigma > 0$ são constantes e $W_t$ representa um processo Browniano.

## Objetivos Teóricos
O trabalho envolve a dedução matemática de diversos métodos de estimação para os parâmetros do modelo:
* **Máxima Verosimilhança (MLE):** Dedução das funções de verosimilhança e log-verosimilhança para obter os estimadores dos parâmetros.
* **Mínimos Quadrados:** Obtenção dos estimadores para $\alpha$ e $\beta$ através de funções estimadoras de mínimos quadrados.
* **Funções Estimadoras de Martingalas:** Dedução de martingalas lineares e lineares aproximadas para a estimativa de parâmetros.
* **Esquema de Euler:** Dedução da função de contraste baseada no esquema de Euler.

## Implementação Computacional
A componente prática do projeto inclui o desenvolvimento de código em **R ou Python** para:
* **Geração de Trajetórias:** Simulação do processo através de densidades de transição e do esquema de Milstein.
* **Simulação Numérica:** Geração de uma trajetória no intervalo $[0, 100]$ com $\Delta = 0.01$, considerando os parâmetros $\alpha=0.1$, $\beta=20$, $\sigma=0.3$ e $x_0=20$.
* **Cálculo de Estimadores:** Implementação de algoritmos para calcular as estimativas dos parâmetros usando a trajetória gerada e otimização numérica para os estimadores de máxima verosimilhança.

---
**Instituição:** NOVA School of Science & Technology - Departamento de Matemática.
