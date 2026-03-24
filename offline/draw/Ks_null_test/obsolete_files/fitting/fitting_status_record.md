# Fit Status Recorder
- plan to only use bin [10, 159]


**version 1:**

 - failed bins: [20,21,22,23,28,29,30,32,33,40,41,42,45,46,48,53,57,62,68,70,75,77,79,83,84,96,98,108,121,134,139,140,152,154,155,157]


```python
 double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":2})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**version 2:**
change the order of chebychev

 - success : [40,42,46,48,57,62,75,83,98,108,140,157]
 - failed bins: [20,21,22,23,28,29,30,32,33,41,45,53,68,70,77,79,84,96,108,121,134,139,152,154,155]


```python
 double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":1})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**version 3:**
change again to 3 order

- success : [23,41,53,70,77,79,84,96,121,]
- failed bins : [20,21,22,28,29,30,32,33,45,68,108,134,139,152,154,155]

- not that kind good : [134,139,152,154,155]
in which the bkg shape forms a bump in signal region, temporarily count as acceptable

```python
 double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":3})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**version 4:**

change order to 1
and the mean range wider -> (0.49761, 0.49461, 0.49961)

- success : [68,108]
- failed : [20,21,22,28,29,30,32,33,45]

```python
double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49461, 0.49961),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":1})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**version 5:**

change order to 3
limit frac range -> (0.7, 0.5, 0.9)

- success : [22,28,33,45] 
- fail : [29,32]
- not so good : [20,21,30]

```python
 double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  "frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":3})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```

**version 6:**

change order to 4
still limit frac

- success : [21,30]
- fail : [20,29,32]

20 use version 5 fit

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  "frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":4})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```


**bin 32:**
not so good
```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":0})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**bin 29:**
not so good
```python
   double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":1})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```
