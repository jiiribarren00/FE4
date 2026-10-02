# Determinación de la relación carga masa del electrón 

## Un poquito de historia

Todavía no se entendían del todo las ondas electromagnéticas, se confiaba en la idea del éter. Ello llevaba a que los rayos catódicos se creyeran "electricidddad en movimiento" una especie de exitaicónd el éter.

### Experimento de Thompson. 
Cátodo, Anodo en un tubo sellado de vidrio. Se exponía el cátodo a una alta diff de potencial y se veía el rayo catódico, hoy sabemos que son electrones. Se ponía una pantalla fluorescente para ver el lugar donde incidía el has de rayos catódicos. 
Se generaba, sobre una parte de la trayectoria del haz, una zona de campos electrico. Bajo la ley de lorenzt $\vec{F}=q(\vec{E}+\vec{v}\times \vec{B})$. Al no tener campo magnético.

## El experimento que haremos nosotros.

> [!CAUTION]
> Para el emisor termiónico, NO superar los 6V en el equipo Pasco o los 10V en el equipo Phywe.

Pasco: $n=130$ y $R=0,15$ m

Phywe: $n=154$ y $R=0,2$ m

Conservación de la energía, entonces $eV=\frac12 m v^2 \Rightarrow v^2= 2\freac{e}{m}V$

COnsiderando la ley de lorentz como una fuerza centrípeta actuando sobre el electrón saliente $F_{centrípeta}=F_{Lorentz}\Rightarrow m \frac{v^2}{r}=evB \Rightarrow \frac{e}{m}=\frac{v}{rB}$.

Considerando la expresión anterior para la velocidad y la forma del campo magnético para las bobinas de Helmholtz: $B=(\frac45)^{\frac32}\mu_0 N i /r$, resulta en $\frac{e}{m}=frac{2V}{B^2r^2}=$

[!NOTE]
Lo que vemos es la desexitación del gas que vemos luego de que los electrones transfieren energía a los electrones del gas. Luego de ello, el electrón seguramente abandona la trayectoria circular.

Tenemos tres variables: $V$, $i$ y $r$. Ordenadamente, mantenagan fija una, varíen una y midan la tercera. Esta parte la repiten para cada una de las tres variables.

# Determinación de $h$

**Objetivo:** determinar la constante de planck $h$.

**Método:** Determinaciones independientes de la irradiancia y temperatura de una boimbilla comercial.

## Parte histórica

## Ley de Planck

## Aproximación de Wien

## Bombilla como cuerpo negro

## Fotodiodo como sensor

El mismo es un diodo convencional, que consta de una interfaz entre dos materiales llamados P y N (positivo y negativo respectivamente), permite el paso de la corriente en un sentido específico y lo impide en sentido contrario. Aunque esto último es una aproximación, en realidad existe una dependencia de la corriente que deja pasar el diodo en función del voltaje.

El diodo tiene un máximo de sensibilidad en $\lambda = 950nm$ de longitud de onda de la luz. Sin embargo, el diodo tiene una curva de sensibilidad con un determinado ancho. Se puede obtener esa curva del Datasheet del componente electrónico. Allí se ve que es una campana (no necesariamente gausiana pero bastante parecida) con una $FWHM = 100nm$.


## Los datos

# Determinación de $\frac{h}{e}$

En este experimento se usa una lámpara de mercurio Hg, pues proporciona varias emisiones en el espectro visible.

## Efecto fotoeléctrico

Pensamos en los fotones como partículas que chocan con un electrón. Entonces, por conservación de la energía:
$$h \nu = \frac12 m v_e^2 + \Phi_0$$
donde $\Phi_0$ es la *función trabajo*, que es el trabajo mínimo necesario para arrancar un electrón del material. Ojo que se confundenla función trabajo, que acá se denominará  $\Phi_0$ pero la notación de la bibliografía puede variar, y el *potencial de contacto*  $\Phi_0 = e \phi_0$, que contiene la misma información pero se le llama potencial en analogía con el potencial eléctrico.

Esos electrones los frenamos con una diferencia de potencial variable $V_f$
$$V_f = \frac{h}{e}\nu-\Phi_0$$

Nos interesará entonces variar las frecuencias y medir $V_f$ para cada $\nu$. Esto nos permitirá calcular $\frac{h}{e}$ y $\Phi_0$.

## Dispositivo experimental

### Lámpara de Mercurio

Empezamos con una lámpara mercurio gaseoso. Esta da una luz blanca pero que en realidad tiene emisión en distintas bandas muy precisas de frecuencia. Podemos ver una simulación y una fotografía de la emisión del Hg en [esta web](https://atomic-spectra.net/spectrum.php?elem=Hg). En [este paper](https://www.researchgate.net/publication/3139924_Radiation_damage_and_light_transmission_studies_on_air_core_light_guides), encontrarán un espectro de emisión muy bonito que les permitirá obtener las longitudes de onda de emisión de una lámpara de mercurio. 

### Red de difracción

Ya conocemos la red de difracción que, por un lado, nos permitirá separar las frecuencias de la luz emitida por el mercurio hacia distintos ángulos.

La red no es perfecta y la difracción es un fenómeno que, en definitiva, es cuántico. Entonces van a tener que ayudar a la red de difracción con un filtro pero pueden probar de hacerlo sin ellos para cuantificar la diferencia.

### Cabezal

El cabezal tiene una pantalla fosforecente que nos permite ver la reflexión de la luz tanto visible como ultravioleta (en realidad esta última se absorbe y reemite en el visible). Esa pantalla, a su vez, tiene una ventana.

Esa ventana permite que la luz pase a una cámara oscura donde se encuentra un bulbo al vacío. Es importantísimo que esté al vacío porque, si no, perderíamos todos los fotoelectrones, siendo estos absorvidos por el aire.

En ese bulbo existe un par de electrodos. Sobre uno de ellos se  

### Circuito amplificador y determinación

Se conecta al par de electrodos 

Probabilidad de las interacciones dobles.

### Filtros

Los filtros no son perfectos. Tener cuidado con eso.

## Papers de caracter histórico

* Paper original de Millikan en donde hace un minucioso trabajo sobre la determinación de la constante de plank https://journals.aps.org/pr/pdf/10.1103/PhysRev.7.355 *Aprecien los diagramas, la dedicación que tienen los graficos y encima hechos a mano*

* Paper de P. Lenard con la primera determinación de $h$ a partir de este método https://zenodo.org/records/1424009

* Paper original de A. Einstein sobre el efecto fotoeléctrico https://onlinelibrary.wiley.com/doi/10.1002/andp.19053220607


# Análisis del Efecto Compton a partir del espectro de emisión del $^{137}$Cs medido con un cristal centellador.

Fundamental para este trabajo la siguiente bibliografía:
[Glenn F. Knoll. Radiation Detection and Measurement. 4.a ed. Wiley, 2010.](https://phyusdb.wordpress.com/wp-content/uploads/2013/03/radiationdetectionandmeasurementbyknoll.pdf)

# Fluorescencia de RayosX

[Paper de Bohr proponiendo las orbitas de energía (momento angular) cuantizado](https://uni-tuebingen.de/fileadmin/Uni_Tuebingen/Fakultaeten/MathePhysik/Institute/IAP/Forschung/MOettel/Geburt_QM/bohr_PhilMag_26_1_1913.pdf)

[Paper de Moseley ajustando rectas a lineas espectrales](https://archive.org/details/londonedinburg6271914lond/page/702/mode/2up)

[Booklet de data espectroscópica](https://amo-csd.lbl.gov/downloads/xraydatabooklet.pdf)
