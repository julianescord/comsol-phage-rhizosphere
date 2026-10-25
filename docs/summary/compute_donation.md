# Escenario de compute donation: Etapa 5 (3D, múltiples beads)

Dimensionamiento técnico del recurso de cómputo que haría falta para ejecutar la Etapa 5, identificada como pendiente en [resumen.html](resumen.html) por falta de hardware adecuado.

## Punto de partida verificable

Etapa 1 (2D axisimétrica, 1 bead): malla de **5,035 nodos**, corre en minutos en la máquina actual (7.5 GB RAM). El propio resumen ya identifica el bloqueo: Etapa 5 no se ejecutó, estimada en "varios millones de grados de libertad", "fuera de alcance" con 7.5 GB, "viable en workstation con ≥32 GB".

## Por qué el salto es tan grande (no es solo "más grande la malla")

1. **Pérdida de simetría de revolución**: 2D axisimétrico → 3D real multiplica nodos ~50-200× para igual resolución.
2. **Capa límite difusiva fina alrededor de cada bead**: la longitud de penetración √(D_suelo/k_inact) ≈ 0.36 mm (Etapa 3) exige elementos sub-décima de mm pegados a cada esfera — no se puede usar malla gruesa lejos y gruesa cerca, hay un gradiente físico real que forzar.
3. **N beads en interacción**: para que el estudio tenga sentido (competencia por *Ralstonia*, heterogeneidad de arreglo espacial, no solo una bead aislada replicada N veces) hace falta resolver varias beads en el mismo dominio, no una sola con condición de simetría.
4. **2 especies acopladas** (fago + *Ralstonia*, como en el modelo combinado) duplican los grados de libertad por nodo.
5. **Transitorio rígido no lineal** (Lotka-Volterra, tolerancia relativa 10⁻⁴-10⁻⁶) con solver directo (MUMPS/PARDISO, lo que usa COMSOL por defecto en dominios acoplados): el costo de memoria de un solver directo en 3D escala ~O(N^1.5), no linealmente. Con N ~ 3-8 millones de DOF eso son fácilmente **80-200 GB de RAM** solo para el factor LU, más margen para el solver iterativo/no lineal.

## Escenario concreto de compute donation

**Qué se necesita:**

- 1 nodo de cómputo con **256-512 GB de RAM** (margen sobre el escalado O(N^1.5) del solver directo, y espacio para 2-3 corridas simultáneas del barrido, no solo una)
- 16-32 cores físicos (COMSOL paraleliza el ensamblaje y el solver iterativo; el directo se beneficia menos pero el pre/post-proceso sí)
- **GPU no es crítico**: COMSOL 6.4 no acelera por GPU el solver de *Transport of Diluted Species in Porous Media*; pedir GPU sería sobre-especificar el escenario
- Licencia COMSOL con módulo **Batch + Cluster Computing** (o LiveLink) para correr sin GUI en el nodo remoto — tan bloqueante como el hardware: sin ese módulo de licencia no se puede ni lanzar el job
- ~200-500 GB de almacenamiento de escritura rápida (SSD) para los archivos `.mph` de cada corrida del barrido (una malla de millones de nodos con resultados transitorios guardados en varios timesteps pesa varios GB por corrida)

**Duración estimada del uso:**

- Una corrida sola (1 punto del barrido: 1 combinación de r_bead, arreglo espacial, Pe) probablemente entre 2-12 h, según cuántas beads y qué tan rígido resulte el transitorio.
- El barrido mínimo útil (análogo en espíritu al barrido 2D ya hecho, sobre r_bead × N_beads × arreglo) son decenas de corridas → **días de cómputo dedicado**, no una corrida única.
- Encaja mejor como una asignación tipo "credits" de un programa académico (semanas de acceso) que como una donación de hardware físico puntual.

## Lo que NO hace falta pedir

- **No** un cluster HPC multi-nodo con MPI: COMSOL en este tipo de problema (un dominio, un solver directo/iterativo dentro de un solo estudio) no escala bien a múltiples nodos, salvo que se reestructure el barrido como jobs independientes por nodo — lo cual sí es razonable, porque cada punto del barrido es *embarrassingly parallel* entre sí, aunque cada corrida individual sigue siendo de un solo nodo.
- **No** GPU.
- **No** almacenamiento masivo tipo petabyte — esto es un problema de RAM y licencia, no de dataset.

## Formulación mínima para pedir el recurso

> Se necesita acceso temporal (2-4 semanas) a 1-2 nodos con 256+ GB de RAM y 16-32 cores, más una licencia COMSOL con Cluster/Batch habilitada, para correr un barrido de ~20-50 simulaciones 3D transitorias con múltiples beads (3-8 millones de DOF por corrida), como validación mínima de la Etapa 5 identificada como pendiente en el proyecto.

Esto es sustancialmente más chico que "un cluster GPU": el cuello de botella real de este modelo es **memoria de un solver directo en 3D + licencia de cómputo por lotes**, no paralelismo masivo ni aceleración por GPU.
