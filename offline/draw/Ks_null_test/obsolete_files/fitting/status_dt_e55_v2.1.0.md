# Fit Status Recorder
- plan to only use bin [10, 149]

**version 1:**

failed_bins = [i+10 for i in range(10)]  + [20,21,22,28,29,72,78,80,84,89,103,104,114,120,123,125,127,130,133,135,137,138,139,142,143,]

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


**version 2:**

- change sigma range
- Ks mass range

success = [11,12,15,18,19,72,78,80,84,89,103,104,125,127,130,133,135,137,139,142,143,]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.520, 100)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.0008, 0.004),  # narror gaussian
                                                  "sigma2":(0.004, 0.0025, 0.006),   # wide
                                                  #"frac":(0.7, 0.5, 0.9)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":1})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```
v2.1 : 2 order bkg
success  = [21,22,114,]
v2.2 : 1 order bkg, Variable("Ks_M", 0.472, 0.525, 106)
success = [120,123,138]
v2.3 : 1 order bkg, Variable("Ks_M", 0.47, 0.515, 90)
success = [13,14,16]
v2.4 : 2 order bkg, Variable("Ks_M", 0.47, 0.515, 90)
success = [10,17]



**version 3:**

- guassian not same mean

- success = [28,29,]

```python
double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                    [PDFSpec("sig", "Ks_M", "double_gauss", 
                                             {"same_mean":False,
                                              "mean1": (0.49761, 0.49561, 0.49861),
                                              "mean2": (0.49761, 0.49461, 0.49961),
                                              "sigma1":(0.002, 0.0008, 0.004),  # narror gaussian
                                              "sigma2":(0.004, 0.0025, 0.006),   # wide
                                              #"frac":(0.7, 0.5, 0.9)
                                             }),
                                    PDFSpec("bkg", "Ks_M", "chebychev",
                                            {"order":2})],
                                    model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```
v5.1 : change to order 1
success : [20]



# Re-Fit
fail = []

**version 1:**
success = []

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":False,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.0008, 0.004),  # narror gaussian
                                                  "sigma2":(0.004, 0.0025, 0.006),   # wide
                                                  "frac":(0.7, 0.5, 1.)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":2})],
                                        model = "nsig[10000,0,2000000]*sig + nbkg[5000,0,200000]* bkg")    
```


