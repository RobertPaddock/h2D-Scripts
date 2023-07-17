# h2D-Scripts

Here we have scripts used to write and analyse h2D simulations. These are based roughly on my hyades-planar-scripts, and I suggest becoming familiar with those first. The 2D versions are not extensively commented, and I have not checked their operation before uploading.

To run, first run the radial and axial meshing scripts (to mesh in the two different dimensions), and then run the file generator to produce the input deck.

Then use the h2DRays script to produce the laser rays. Take the output from this, and copy it into the appropriate point in the input deck.

It can then be analysed using the h2D Reader file.

These scripts likely will require some work, but may be a useful start point for anyone hoping to do planar sims in h2D.