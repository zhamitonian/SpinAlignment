# Hadronic Selection

## Main flow

### Paricle lists :

**track** : 
- $p_t > 0.1$, $dr < 2$, $dz<4$. 
- TrkList, Mdst_charged_Manager. 

**gamma** : 
- unmatch to any track, $E>0.1$. 
- GamList, Mdst_gamma_Manager. 

**cluster** :
- unmatch, $E>0.1$, $17^{\circ} <\theta < 150^{\circ}$
- EclListCMSel, Mdst_ecl_Manager. 

### Variables :

**HadronA:**

- $\sqrt{s}$ = 
- $N_{trk} \geq 3$ ,(length of TrkList)
- C.M. variables, boost with $\pi$ and $\gamma$ mass hypothesis
- $E_{vis}^{\star} \geq 0.2*\sqrt{s}$ , ($\sum_{GamList} E^{\star}$ + $\sum_{TrkList}E^{\star} $ )
- $|\sum p_z^{\star}| \leq 0.5*\sqrt{s}$, ($\sum_{GamList,TrkList} p_z^{\star}$ ) 
- $ 0.025 \leq E_{sum}/\sqrt{s} \leq 0.9$, ($E_{sum} = \sum_{EclListCMSel} E^{\star}$)
- $PrimeR \leq 1.5, |PrimeZ| \leq 3.5$

**new HadronA:**
- $ 0.1<E_{sum}/\sqrt{s} < 0.8$
- $Necl_{barrel} > 1$ , (number of EclListCMSel with $-0.7 <\cos\theta < 0.9$ )

evtinfo func 中的值是用的之前的这些list以及上面描述的变量定义。

**HadronB:**

