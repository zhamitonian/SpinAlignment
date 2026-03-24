# Fit Status Recorder
- [10, 159]

**version 1:**

fail : [10,11,12,13,14,15,16,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,52,54,66,72,79,86,91,101,124,141,142,145,150,151,152,154,157,]

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

**version 2:**

apply 2nd order chebychev

fail : [11,13,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,52,54,66,72,79,86,101,141,142,145,150,151,152,154,157,]

success : [10,12,14,15,16,40,124,]

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

**version 3:**

apply 3rd order chebychev

fail : [13,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,44,46,47,48,141,142,145,150,151,152,154,157,]

success : [11,41,42,43,45,52,54,66,72,79,86,101,]

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

change x range (to avoid bkg get 0)

fail : [,20,24,25,26,27,29,32,33,34,35,36,37,44,46,47,48,141,142,145,150,151,152,154,157,]
success : [13,18,19,21,22,38,]
not so good : [25,26,32,33,34,35,46,48,]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
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

**version 5:**

2 order
fail : [20,24,25,26,29,32,34,35,36,37,141,142,145,150,151,157,]
success : [27,33,44,46,47,48,152,154,]
nos so good : [20,29,]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
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

**version 5:**

some mean

fail : [24,25,26,29,32,34,35,36,141,142,145,150,151,157,]
success : [20,37]
not so good but acceptable : [20,29,]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
    #double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":False,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "mean2": (0.49761, 0.49461, 0.49961),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":3})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```

**extra1:**
[25,26,32,34,35] use version 4
[29] use version 5

**version 6:**

polynomial and limit frac

no so good : [141,142,145,150,151,157] but use this
fail : [24,36]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  "frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "polynomial",
                                                {"order":2})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")   
```

**version 7:**

not good : [24,36] but use it

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  "frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "polynomial",
                                                {"order":3})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")   

```