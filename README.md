# Wavelet Transform for dummies(GUI)

Wavelet Transform for dummies is a MATLAB based GUI aimed at helping people with little to no programming knowledge perform analysis using Wavelets

This program can read a signal, calculate it's wavelet transform and plot it's average amplitude/power. One can also save the plots, Wavelet Transforms and average vectors.

This release currently needs to be run on MATLAB. An installer which does not need MATLAB will be put up soon.

## Getting Started
Download all the files in this repository along with all the ones in the auxilaryFunctions repository. Add both to path in MATLAB and simply run Copied.m

## Instructions
#### Importing a signal
Click on File to load your signal in csv or .mat format

#### Selecting Data
- Use the zoom tool to select your data. Once you have chosen it click on Refresh
- If you would like to type in the data limits. Simply type it into the box labelled Xlim

#### Wavelet Transform Options
- Sampling Frequency: Rate at which the signal was sampled

- Max Frequency: Maximal frequency for which to calculate Wavelet Transform. Default value is sampling_frequency/2

- Min Frequency: Minimal frequency for which to calculate Wavelet Transform. Default value is the minimal frequency for which at least one WT coefficient is determined

- Central Frequency: Central frequency of the wavelet

- Under Sampling Rate: The number indicates every nth point that is used for plotting. WARNING: Going too high can lead to erroneous results

- Amplitude/Power: Plots the amplitude or the power plot of the wavelet transform

#### Advanced Options
If you don't know what you're doing. DO NOT TOUCH THIS
- Wavelet Type: The type of wavelet used for the wavelet transform

- Preprocess: Subtract the 3rd order polynomial fit and then bandpass the signal in the band of interest [fmin,fmax]

- CutEdges: Determine if WT coefficients should be set to NaNs out of the cone of influence

#### Buttons
- Wavelet Transform: Calculates and plots the wavelet transform

- Plot: Used to refresh the plot without a need to need to calculate the wavelet transform again

- 2D Plot: Rotates the view point of the surf plot to the xy-plane 

- Intervals: Plots a dashed(black) line on the plots to facilitate easier viewing 

## What is a wavelet transform?
Okay, so you obviously you know what a fourier transform is right? Well if you don't, in simple terms

> ### Give a Fourier Transform a Krabby Patty and it'll give you the list of ingredients

(Too bad Plankton didn't know this)

Pretty cool huh? Now what if I told you the Wavelet Transform is **_better_** than a Fourier Transform?

#### Why is it better?
> The main difference is that wavelets are localized in both time and frequency whereas the Fourier transform is only localized in frequency

If that made no sense to you. You've probably never heard of the Uncertainity principle.

---

#### The Heisenberg Uncertainty Principle
So basically it turns out that the **position** and **momentum** cannot be measured **simultaneously** with an **arbitrarily high** degree of precision

More information over here: [Uncertainity Principle][1]

---

##### Okay so what?

You see, in the fourier transform we go to the limiting case of the uncertainity principle viz. our bandwidth is 0. Now this means that we know the size of the frequency component **exactly** but we have no information about the temporal spread of the spectrum.

Now the wavelet transform is smart. It strikes a balance between the two in order to retrieve both temporal and spatial information, neither of them exact of course, but sometimes that's all you need.

#### When to use a wavelet tansform

Like the saying goes. 
> You don't need a blow torch to light a spliff

An example should clear things up and help you decide what to use 

## Decrypting the Krabby Patty

So let's say Plankton is trying to steal the secret formula(again) 

Now Plankton snuck into the kitchen and planted a detector in the kitchen. This detector basically reads whatever is there on the stove and records it in the form of an digital signal.

Now that he knows about Fourier Transforms he can put his diabolical plan into action. He pulls out his FT-ingredient analysis-inator err... let's just call it a cool device. This device calculates the fourier tranform and displays what ingredients were used.

Fast forward to Plankton with readings in hand. Plankton is elated and goes off on his monologue about how he's going to run Mr.Krabs out of business. Plankton rushes to the grocery store and gets all the ingredients. It's not until he comes back that he realises he's messed up pretty bad. He knows what's **_in_** a Krabby Patty but doesn't know how to put it together or even in what order!

Sad and depressed, Plankton spends his days in depression looking at pictures of Krabby Patties  on the internet. 

When one day he comes across this book. [The Illustrated Wavelet Transform Handbook][2]

Now Plankton is armed with the wavelet transform. He goes through with the same plan except this time uses his WT-ingredient analyser.

**insert diabolical laughter**

He interprets the signal like this 

- At the first parameter set w let's say t(w) is **roughly** a minute or so since Spongebob started cooking and k(w) is **something similar to bread** like a bun. This tells Plankton that up until a minute or so he just needs to heat the bun.
- He just repeats this for the rest of the signal and voila he has an **approximate** but hopefully close enough rip off of the Krabby Patty.

Just raise an issue if you have any questions or suggestions

e-mail: guruvamsi.policharla@gmail.com

#### Credits
Thanks to everyone in the Nonlinear Biomedical Physics Group, Lancaster University for their constant support!

[1]: http://hyperphysics.phy-astr.gsu.edu/hbase/uncer.html       "Uncertainity Principle"
[2]: https://books.google.co.uk/books?id=RUSjIMQACQQC&source=gbs_similarbooks "Addison_Book"

