# seph
Asteroid ephemeris demo

Someone wanted an ephemeris program and said I should be able to do it
because it's all in Meeus.  So, here is a little something I hacked up.
It uses packages on this account from the meeus repo as well as the
sexagesimal repo recently split off from meeus and the mpcformat repo,
for reading the file MPCORB.DAT.

It did take a little hacking.  In mpcformat, I found unfinished code
to parse the MPC packed epoch format, so I patched that up and pushed it.

Then it turned out that while the meeus packages were sufficient for
computing RA and dec, there were parts missing for computing apparent
magnitude.  The math is all in Meeus's book, but this is one of the cases
where there was no worked example in the book and so I had not implemented
functions in the package.  (That was a criterion for inclusion.  With no
worked example, there was no test data provided and so no simple correctness
test.  I wanted a test for each included function.)

These parts missing from meeus are in the seph.go source file here.
There is a small function for computing magnitude but it takes some
quantities as inputs that are not computed or not exported by other
meeus functions.  For those, I cut and pasted a couple of functions
from meeus and modified them.

The initial commit demonstrates a capability, that's all.
