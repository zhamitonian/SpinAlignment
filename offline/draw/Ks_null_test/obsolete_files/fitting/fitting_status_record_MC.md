# Fit Status Recorder
- [10, 159]

**version 1:**
two gussian with not same mean value

fail : [11,12,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,69,76,81,82,87,92,94,98,101,104,105,106,107,108,141,143,144,145,148,149,151,153,156,159]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":False,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "mean2": (0.49761, 0.49461, 0.49961),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":1})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    

```

**version 2:**
change to 2 order bkg

fail = [11,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,82,92,143,144,149,151]
success = [12,69,76,81,87,94,98,101,104,105,106,107,108,141,145,148,153,156,159]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":False,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "mean2": (0.49761, 0.49461, 0.49961),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":2})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```

**version 3:**
use same mean

fail = [11,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,82,144,151]
success = [,40,92,143,149]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 100)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  #"mean2": (0.49761, 0.49461, 0.49961),
                                                  "sigma1":(0.002, 0.001, 0.01),  # narror gaussian
                                                  "sigma2":(0.01, 0.0001, 0.02),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":2})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```

**version 4:**

