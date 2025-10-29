# fits1

This exercise will begin to introduce non linear data fitting.  The examples will be based on the ROOT libraries.  Similar functionality for everything we see can be found in the numpy + scipy + lmfit + matplotlib modules. But the examples we will follow over the next few projects offer several advantages:
* far more natural histogramming tools
* completely automated fitting options ("one liners")
* convenient methods for defining custom functions and controlling the convergence of the fit
* detailed and consistent documentation
* and a low level interface to modify the objective function, running in fully optimized compiled code

You are welcome to modify the provided code for your projects and to use other packages.  Where applicable alternate examples will be included. 

* **fit1.C**: C++ function to generate random data according to a normal distribution with mean=20, sigma=10. <br> A fit is performed using the built in Gaussian model in ROOT.  Then the parameter values, their uncertainteis, and the p-value for the fit are extracted.  To run this code type ```root fit1.C``` or if you are already running ROOT, type ```.X fit1.C```  
* **fit1.py**: The same code using the python interface, run this example using ```python fit1.C```.
* For a contrast see **fit1mpl.py** for a version using matplotlib+scipy.  
* readhist.C(py):  Examples for reading the histogram files given in this example 
* ParamUnceratinties.ipynb : a guided tutorial towards most of what you will be coding in this week's exercise.
* LLexample.ipynb : a notebook giving an example for calculating (N)LLs
* TH1hist2Numpy.ipynb : an example for converting a ROOT histogram to numpy arrays

Note that from ROOT you can type ```new TBrowser()``` or in Python r.TBrowser() to get a graphical browser that allows you to look at what's contained in the TFiles.


Exercise 1.
**How does the 1 sigma width of this distribution compare to the typical size of the uncertainty reported for this fit parameter?**
**Does this seem reasonable? Discuss.  (eg. add a comment to you Readme.md)**
The sigma width for the chi square., 0.2125, and the error mean is 0.3257. I think it's pretty reasonable, considering how small the variation in chi2 is.

Exercise 2.
**How do your results compare to the expected values in each case? How do the distributions of the parameter values from the fits compare to the estimated uncertainty on the fit parameters?**

The NLL models decay much faster and go to zero unlike the chi2 away from the fit parameters. Chi2 has many more non-zero bins further way from the center
answer for the fit.

Exercise 3.
**Estimate the "p-value" of your original fit an give this value in your Readme.md.**
The p-value I got 0.458. 

Exercise 4.
Generally the NLL seems to produce narrow contours compared to the chi2 value.