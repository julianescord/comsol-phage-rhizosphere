<!-- doc2md origen="/media/julianescord/DATA/Documentos/Antigravity/comsol/docs/summary/resumen.html" sha256="63b9b62d8a5eedc8" -->

# Liberación controlada de bacteriófagos en la rizosfera

Metodología y hallazgos de un modelo de elementos finitos que adapta la técnica de Dosmar et al. (2021) -desarrollada para fármacos oculares de liberación sostenida- a la entrega de bacteriófagos anti- *Ralstonia solanacearum* desde una micro-esfera de alginato hacia la rizosfera.

**Software** COMSOL Multiphysics 6.4 **Construcción** programática, Python + `mph` **Geometría** axisimétrica 2D **Etapas** 6, todas verificadas **Actualizado** 24 jul 2026

6

etapas construidas y verificadas de forma independiente

0.36 mm

longitud de penetración difusiva del fago en suelo

R₀≈4.3

número reproductivo básico del fago sobre *Ralstonia*

300 µm

vehículo óptimo una vez incluida la replicación del fago

&lt;1-2%

error de cierre del balance de masa en cada etapa

1. [Alcance Objetivo y enfoque](#objetivo)
2. [Método Metodología](#metodologia)
3. [Etapa 1 Difusión pura](#etapa1)
4. [Etapa 2 Suelo real + inactivación](#etapa2)
5. [Etapa 3 Captación radicular](#etapa3)
6. [Etapa R Replicación del fago](#etapaR)
7. [Etapa 4 Flujo de agua](#etapa4)
8. [Combinado Flujo + replicación](#combinado)
9. [Síntesis Discusión consolidada](#discusion)
10. [Honestidad Límites del modelo](#limitaciones)
11. [Cierre Estado y reproducibilidad](#estado)

La numeración refleja el orden real de construcción, incluida una reordenación deliberada (la etapa de raíz se adelantó a la de replicación).

· Alcance

## Objetivo y enfoque

Predecir cuánto bacteriófago activo llega desde un vehículo de alginato hasta la superficie radicular, y qué combinación de tamaño de vehículo y velocidad de liberación maximiza la supresión de *Ralstonia solanacearum* .

El fago no se comporta como el fármaco pasivo del paper de referencia: además de difundirse y perder actividad en tránsito, se **amplifica** al infectar a su hospedador. Esa diferencia -inactivación como sumidero, replicación como fuente- es el eje de todo el trabajo, y resultó decisiva: la respuesta a "¿qué tamaño de vehículo conviene?" se **invierte por completo** según qué procesos biológicos se incluyan.

El proyecto se construyó como una secuencia de seis modelos de complejidad creciente, cada uno verificado antes de pasar al siguiente. Ningún resultado se reporta sin una comprobación independiente -solución analítica, invariante de conservación, caso límite o comparación contra una ODE de referencia.

· Método

## Metodología

Geometría axisimétrica 2D en *Transport of Diluted Species in Porous Media* ; construcción y verificación 100% programáticas para que cada modelo sea reproducible desde cero.

La rizosfera se representa como un dominio de suelo con una micro-esfera de alginato embebida (revolucionada en torno al eje de simetría). Todos los modelos se construyen con un paquete Python ( `mph` ) que controla la API Java de COMSOL 6.4 directamente -geometría, física, malla, estudio y post-proceso quedan fijados en el script, no en clics de la interfaz gráfica-, de modo que cualquier resultado se puede reconstruir ejecutando un archivo.

### Variables de entrada (parámetros del modelo)

|Símbolo|Descripción|Valor / rango|Unidad|
|---|---|---|---|
|r_bead|Radio de la micro-esfera de alginato solo escala micro|25 - 1000 (barrido)|µm|
|L_root|Distancia bead → superficie radicular|0.1 - 5|mm|
|D_water|Difusividad del fago en agua libre (Stokes-Einstein, partícula ~200 nm)|2.2×10⁻¹²|m²/s|
|D_bead|Difusividad aparente del fago dentro del gel de alginato|2.2×10⁻¹² - 10⁻¹⁵ (barrido)|m²/s|
|D_host|Difusividad de *Ralstonia* (sésil, biofilm)|1×10⁻¹⁶|m²/s|
|ε\_soil / ε\_bead|Porosidad del suelo / fracción de agua del gel|0.45 / 0.98|-|
|k_inact|Tasa de inactivación del fago|0.5|1/d|
|c₀|Concentración inicial de fago en la bead|1|mol/m³|
|Vmax\_root , Km\_root|Captación radicular saturable (Michaelis-Menten)|1×10⁻⁹ mol/(m²·s), 5×10⁻³ mol/m³|-|
|b_burst|Fagos liberados por cada lisis (burst size)|50|-|
|k_inf|Tasa de infección fago- *Ralstonia*|5×10⁻⁷|m³/(mol·s)|
|r\_host , K\_cap , H₀|Crecimiento logístico, capacidad de carga y densidad inicial de *Ralstonia*|0.5 1/d, 1 mol/m³, 1 mol/m³|-|
|Pe_ref|Número de Péclet del flujo de agua de fondo|0 - 100 (barrido)|-|
|t_end|Horizonte de simulación|14 - 21|d|

Todos los barridos de r\_bead ejecutados en COMSOL están en escala **micrométrica** (25-1000 µm). La escala nanométrica se evaluó por separado, solo de forma analítica - ver el recuadro siguiente.

Caso límite · por qué no se simuló la escala nanométrica

La Etapa 1 ya había mostrado que, con difusividad uniforme, el problema es **auto-similar** : el tiempo característico de liberación escala como τ = a²/D , así que basta este cálculo cerrado -sin necesidad de FEM- para saber si una nanoesfera es viable antes de gastar una corrida en ella.

|Radio|τ = a²/D (D = 2.2×10⁻¹² m²/s)|Régimen|
|---|---|---|
|50 nm|1.1 ms|liberación instantánea|
|100 nm|4.5 ms|liberación instantánea|
|200 nm|18 ms|liberación instantánea|
|500 nm|114 ms|liberación instantánea|
|1000 nm (1 µm)|455 ms|liberación instantánea|

En todo el rango 50 nm - 1 µm, τ &lt; 0.5 s : la carga sale de la esfera antes de que el experimentador termine de pipetear, muy lejos del horizonte de 14-21 días que exige el diseño. Y aunque el tiempo fuera aceptable, hay una restricción física previa que ya descarta la mayoría de este rango: el vehículo debe ser mayor que su carga. Los fagos jumbo anti- *Ralstonia* miden **~250-290 nm** de longitud total, así que no caben en una nanoesfera de 100 nm; solo un podovirus pequeño (~70 nm) encajaría, y muy al límite. Por eso el barrido en COMSOL arranca en 25 µm y no explora nanómetros: la escala nano queda excluida por dos argumentos independientes -geométrico y cinético- antes de llegar a necesitar una simulación.

### Variables de salida (campos y métricas derivadas)

|Variable|Significado|Uso en el análisis|
|---|---|---|
|c / cP|Concentración de fago en el fluido de poro|Campo primario; base de todos los balances de masa|
|cH|Concentración de *Ralstonia*|Dinámica del patógeno (etapas con replicación)|
|Mcapt / Mout|Fago capturado o salido por la raíz (ODE global, integrado por el solver)|Balance de masa; fracción entregada|
|Minact|Fago inactivado acumulado (ODE global)|Balance de masa; pérdidas en tránsito|
|Mprod / Hinf|Fago producido por lisis / *Ralstonia* consumida (ODE global)|Invariante de acoplamiento exacto: Mprod = b\_burst · Hinf|
|carga acumulada|(1/t_end)·∫ H(t)/H₀ dt|Métrica **robusta** de biocontrol (insensible a la fase del ciclo)|
|supervivencia final|H(t_end)/H₀|Métrica **frágil** : solo válida como referencia cruzada|
|amplificación|pico(cP)/c₀|Cuánto se multiplica el fago tras la replicación|
|Péclet Pe|u·L / D_e,suelo|Régimen de transporte: difusivo (Pe≪1) vs. advectivo (Pe≫1)|
|residual de balance|Error de conservación de masa en cada instante|Verificación numérica universal, en todas las etapas|

### Estrategia de verificación

Cuatro herramientas, aplicadas según lo que cada etapa permite comprobar: **(1)** solución analítica de Crank para una esfera con concentración inicial uniforme, autoverificada por conservación de masa antes de usarse como referencia; **(2)** configuraciones degeneradas con solución cerrada (porosidad=1 y tortuosidad=1 reducen el transporte en medio poroso a la 2.ª ley de Fick pura); **(3)** comparación contra una ODE de Lotka-Volterra integrada de forma independiente ( `scipy.solve_ivp` ); **(4)** invariantes de conservación exactos (p. ej. cada lisis produce exactamente b\_burst fagos nuevos) y balances de masa de varios términos, calculados con acumuladores GlobalEquations integrados por el propio solver.

1 Etapa 1

## Difusión pura en suelo homogéneo

Configuración deliberadamente reducida a la 2.ª ley de Fick, para poder validarla contra una solución analítica exacta.

Con porosidad = 1 y el modelo de tortuosidad en τ\_F = 1, la interfaz de medios porosos se reduce exactamente a difusión fickiana (D\_e = D), lo que permite comparar el campo FEM completo -en los 5035 nodos de la malla, no solo en una línea- contra la solución de Crank para una esfera con carga inicial uniforme liberando en un medio infinito.

Validación FEM vs. analítico. Izquierda: perfiles radiales en 4 instantes (puntos = FEM en todos los nodos, línea = Crank). Derecha: error máximo y RMS normalizados frente al umbral de aceptación (2%).

<!-- image -->

|Comprobación|Resultado|Veredicto|
|---|---|---|
|Error máximo FEM vs. Crank, todo el transitorio|4.2×10⁻³·c₀|Pasa|
|Error RMS a t = 24 h|1.2×10⁻⁴·c₀|Pasa|
|Conservación de masa en el tiempo|exacta|Pasa|

El offset geométrico de +0.13% en la masa es discretización del contorno curvo, no un error de física.

2 Etapa 2

## Difusividad propia del gel, suelo real e inactivación

Se activan tres elementos que el paper de referencia no necesitaba: D\_bead distinta de D\_water , porosidad real del suelo con difusividad efectiva de Millington-Quirk, e inactivación del fago como sumidero de primer orden.

No se incluyó sorción de equilibrio (K\_d lineal): un fago de ~200 nm es un **coloide** , no un soluto disuelto, y su retención en suelo es adhesión/filtración cinética -casi irreversible- y no un reparto instantáneo. Queda modelada como un sumidero cinético adicional, desactivado por defecto.

### Verificación por partes

|Comprobación|Configuración de prueba|Error|Veredicto|
|---|---|---|---|
|Inactivación|Medio uniforme, τ=1, contra Crank·exp(−k·t)|4.2×10⁻³·c₀|Pasa|
|Millington-Quirk|Medio uniforme, ajuste libre de D efectiva|4.3×10⁻³·c₀|Pasa|

El ajuste libre del segundo check recuperó un exponente de **0.342** frente al 1/3 teórico de Millington-Quirk saturado -una confirmación empírica, no asumida, de qué fórmula aplica realmente el solver.

Barrido r\_bead × D\_bead. Máximo de fago activo entregado al suelo (%). El resultado es monótono : la esquina superior izquierda (bead pequeña, gel abierto) siempre gana.

<!-- image -->

Hallazgo

Sin ningún sumidero que premie la sostenibilidad, la entrega de fago activo es monótona: **gana siempre el vehículo más pequeño y más difusivo** . No hay ningún argumento, en este modelo, para una liberación sostenida.

3 Etapa 3

## Captación radicular saturable

Sumidero de Michaelis-Menten en la frontera radicular, a distancia L\_root de la bead -el proceso con más posibilidades de romper la monotonía de la Etapa 2.

Esta geometría no tiene solución cerrada, así que se verificó con un invariante exacto: en todo instante, M₀ = activo(t) + captado(t) + inactivado(t) . El balance destapó un detalle real de COMSOL -el nodo FluxBoundary trae species = 0 por defecto, y sin activarlo la condición de frontera **no se aplica** -; se detectó porque el balance daba lo físicamente imposible (todo inactivado y, además, masa capturada).

|Comprobación|Resultado|Veredicto|
|---|---|---|
|Balance de masa de 4 términos|2.1×10⁻⁵|Pasa|
|Límite V_max→0 (debe reproducir la Etapa 2)|captado 0.000%|Pasa|

Escala física

La longitud de penetración √(D\_suelo/k\_inact) ≈ 0.36 mm predice el colapso de la captación con la distancia: 51% a 0.3 mm → 0.09% a 5 mm. **La bead debe estar sub-milimétricamente cerca de la raíz** para que el fago llegue vivo.

Barrido con captación saturable (L = 0.5 mm). El máximo sigue en la esquina (86.8%, bead pequeña y difusiva): la saturación no basta para romper la monotonía.

<!-- image -->

R Etapa de replicación

## *Infección fago- Ralstonia : el aporte de novedad*

El fago deja de ser pasivo: se amplifica al infectar a su hospedador. Sistema depredador-presa espacial de dos especies ( cP fago, cH *Ralstonia* ), cinética Lotka-Volterra.

*Ralstonia* coloniza una banda rizosférica cerca de la raíz, en su capacidad de carga; el fago liberado desde la bead debe alcanzarla y amplificarse. Con los parámetros actuales, el número reproductivo básico es R₀ ≈ 4.3 .

|Comprobación|Referencia|Error|Veredicto|
|---|---|---|---|
|Caso bien mezclado vs. ODE|Lotka-Volterra integrada con `scipy`|1×10⁻³|Pasa|
|Invariante de acoplamiento|Mprod = b_burst · Hinf|1.5×10⁻¹⁵|Pasa|
|Límite H₀→0 (sin hospedador)|decaimiento exp(−k_inact·t)|1×10⁻³|Pasa|

El sistema amplifica ~170× y es numéricamente rígido ( *stiff* ): la tolerancia por defecto del solver bastaba para un ~1% de error en el pico; se corrigió apretando la tolerancia relativa a 10⁻⁴-10⁻⁶.

### El hallazgo central del proyecto

El sistema es un **ciclo depredador-presa con retardo de transporte** : *Ralstonia* crece, el fago la aplasta, el patógeno **rebrota** . Por eso la supervivencia medida en un instante fijo es una métrica frágil -depende de en qué fase del ciclo se mire-; la métrica robusta es la **carga acumulada de patógeno** , ∫H(t)/H₀ dt en toda la ventana de simulación.

Barrido con replicación. Izquierda (métrica robusta): la carga acumulada mejora con la bead grande y de gel lento . Derecha (métrica frágil): la supervivencia a t=21 d muestra un óptimo interior distinto -evidencia de que la fase del ciclo importa.

<!-- image -->

Inversión de la conclusión

Una bead de 25 µm amplifica el fago **46 000×** pero controla peor (65% de supervivencia) que una de 75 µm (1432×, 44%). Más amplificación no es mejor biocontrol: lo que importa es **sostener la presión** sobre el rebrote, no el pico. Dos de tres métricas -incluida la robusta- favorecen beads **grandes y de liberación lenta** : la conclusión de las etapas 2-3 se invierte.

4 Etapa 4

## Flujo de agua y transporte advectivo

¿El flujo de agua del suelo rescata al fago de la inactivación en tránsito hacia la raíz?

Salvedad metodológica

El plan original -transpiración convergente resuelta con la ley de Darcy- no se pudo modelar limpio: la "raíz esférica envolvente" heredada de las etapas 2-3 solo admite un flujo convergente incompresible con una fuente singular en el centro, y regularizar esa singularidad rompe el balance de masa (verificado: el campo convergente da errores de ~1000×). Se optó por un **flujo de fondo uniforme** (percolación que atraviesa la rizosfera), que sí conserva masa exactamente. Una transpiración convergente real exigiría rediseñar la geometría.

|Comprobación|Resultado|Veredicto|
|---|---|---|
|Balance de masa (varios niveles de Pe)|≤ 1.1×10⁻²|Pasa|
|Límite Pe→0 (difusión pura)|1% entregado|Pasa|
|Entrega monótona creciente con Pe|-|Pasa|

Efecto del número de Péclet. Sin flujo, el 100% del fago se inactiva antes de llegar; con un flujo modesto (Pe≈5) el 80% llega vivo. La transición ocurre en Pe≈2.5.

<!-- image -->

### El flujo rescata al fago de la inactivación en tránsito

Fracción de fago entregada viva a la raíz: **Pe=0 → 1%** (todo se inactiva) · **Pe=5 → 80%** · **Pe=20 → 96%** · **Pe=100 → 99.7%** . Esto relaja la restricción de la Etapa 3: la bead ya no necesita estar pegada a la raíz si hay flujo que la transporte.

+ Modelo combinado

## Flujo + replicación

¿El flujo de agua cambia la conclusión sobre liberación sostenida que dio la Etapa R?

El fago se advecta; *Ralstonia* , sésil, no debe hacerlo -con Pe≈20 el flujo la desplazaría decenas de milímetros en 21 días, algo sin sentido biológico-. Como en la interfaz de transporte la velocidad es común a todas sus especies, se usaron **dos interfaces acopladas** : una para el fago (con advección) y otra para *Ralstonia* (sésil), unidas por las reacciones cruzadas de Lotka-Volterra.

|Comprobación|Resultado|Veredicto|
|---|---|---|
|Invariante Mprod = b·Hinf (con advección)|1.5×10⁻¹⁵|Pasa|
|Límite Pe→0 reproduce la Etapa R|85.7% vs. 86.9%|Pasa|
|Balance de masa del fago|1.25×10⁻²|Pasa|

Barrido Pe × r\_bead. La bead de 300 µm (recuadros) gana en las tres filas de flujo: el flujo refuerza , no invierte, la preferencia por liberación sostenida.

<!-- image -->

Sinergia bead-grande + flujo

Para la bead grande, más flujo siempre ayuda (75%→61%→57% de carga). Para beads pequeñas hay un óptimo: **Pe=5 supera a Pe=20** , porque un flujo demasiado intenso barre el fago fuera del dominio antes de que amplifique. La liberación sostenida es más robusta al flujo, porque repone continuamente el fago que el agua arrastra.

· Síntesis

## Discusión consolidada

La misma pregunta de diseño -¿qué vehículo conviene?- recibe respuestas opuestas según qué física se incluya.

|Física incluida|Vehículo óptimo|Dirección|
|---|---|---|
|Difusión pura (Et. 1)|-|Problema autosimilar: no discrimina tamaño|
|+ D_bead propia, inactivación (Et. 2)|Pequeño, gel abierto|Liberación rápida|
|+ captación radicular saturable (Et. 3)|Pequeño, junto a la raíz|Liberación rápida (sin cambio)|
|+ replicación del fago (Et. R)|Grande, gel lento|Liberación sostenida - se invierte|
|+ flujo de agua (combinado)|Grande + flujo moderado|Liberación sostenida - se refuerza|

La conclusión final **valida la intuición original del proyecto** -micro-beads de alginato de liberación sostenida- pero con un matiz que solo el modelo podía revelar: esa recomendación depende enteramente de que el fago se replique en la rizosfera (R₀&gt;1). Con transporte puramente pasivo, el modelo habría recomendado justo lo contrario.

! Honestidad

## Límites del modelo

- **La conclusión depende de la métrica temporal.** El sistema depredador-presa tiene rebrote; una lectura en un instante fijo puede invertir la lectura de la carga acumulada. Se reportan ambas siempre que es relevante.
- **Los parámetros de la cinética fago-** ***Ralstonia*** **son placeholders** de literatura genérica, no medidos para este sistema. La *dirección* del hallazgo (sostenida &gt; rápida, con replicación) es robusta a variarlos; la *magnitud* exacta (300 µm, R₀=4.3) no lo es.
- **El flujo de agua es un campo uniforme impuesto** , no una solución de Darcy/Richards completa. La geometría de raíz envolvente -óptima para difusión y reacción- es incompatible con un flujo convergente incompresible sin una fuente singular; modelar transpiración convergente real exige rediseñar la geometría.
- **La Etapa 5 (3D, múltiples beads) no se ejecutó.** Una estimación de costo con el hardware disponible (7.5 GB de RAM) indica varios millones de grados de libertad para geometrías con más de una bead -fuera de alcance en esta máquina; viable en una workstation con ≥32 GB, o como validación mínima de una sola bead.

✓ Cierre

## Estado y reproducibilidad

Seis modelos construidos y verificados, cada uno como un commit independiente con su justificación documentada.

3cadb80 Estado inicial: infraestructura COMSOL + Etapa 1 validada

e754503 Etapa 2: liberación con D propia, suelo poroso e inactivación (verificada)

933371b Etapa raíz: captación radicular saturable (verificada); no crea óptimo interior

622579d Etapa replicación: infección fago-Ralstonia; invierte la conclusión de diseño

44c00d2 Etapa 4: flujo de agua y advección (campo uniforme); el flujo rescata al fago

d48931b Combinado flujo + replicación: el flujo refuerza la liberación sostenida

### Próximos pasos

El paso de mayor impacto es conseguir **parámetros reales de literatura** (burst size, tasa de infección, crecimiento y densidad de *Ralstonia* , difusividad del fago en alginato): los barridos ya están listos para recibirlos y convertir esta dirección cualitativa en una recomendación cuantitativa defendible. En paralelo quedan abiertos un refinamiento de la cinética de infección (clase infectada con latencia explícita) y la extensión a 3D cuando haya hardware adecuado.