# Fit Status Recorder
- plan to only use bin [10, 149]

failed = []

**version 1.0:**
success = [60,70,71,73,74,75,76,77,78,79,80,83,84,85,86,89:146,148,149]

```python
    double_gauss_config = FitDefinition([Variable("Ks_M", 0.47, 0.525, 110)],
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
v1.1 : change x range to (0.47, 0.52, 110)
success += [12,13,14,15,16,17,]
v1.2 : change bkg to 2nd order
success += [50,51,52,64,67,69,81,82,87,88]
v1.3 : change bkg to 3rd order
success += [11,30,39,40,49,57,58,61,63,65,66,72,]
v1.3 : not use same mean , "mean2": (0.49761, 0.49461, 0.49961),
success += []
v1.4 : v1.3 + 2nd bkg
success += [10,]
v1.5 : v1.3 + 3rd bkg
success += [20,25,28,29,31,35,38,53,54,59,]
v1.6 : v1.5 + change x range to (0.47, 0.52, 110)
success += [18,19,]
v1.7 : change x to  (0.475, 0.52, 90) + v1.3
success += [47,48,62,68,]
v1.8 : v1.7 + 2nd bkg
success += [41,]
v1.9 : change x to  (0.475, 0.52, 90)
v1.10 : v1.9 + 3rd bkg
success += [32,42,45,]
v1.11 : v1.3 + v1.9, "frac":(0.7, 0.5, 0.9), 2nd bkg
success += [21,22,23,26,37,44,55,]
v1.12 : v1.3 + frac":(0.7, 0.5, 1), 3rd bkg
success += [24,27,56]
v1.12 : v1.3 + frac":(0.7, 0.5, 1), 2rd bkg
success += [33,34,36,43,46,147]

