Bake runs a set of code a number of times, each time having a different set of input values.

Bake reads a bake parameter file. A bake parameter file is a list of keys which
have one or more values. Then, bake does something for each combination of
values.  It's useful, for doing repetitive find-and-replace operations,
wrangling data out of oodles of really-similar-but-sublty-
different-in-two-variables-but-not-the-other-five sets of data, doing the
accounting on submitting jobs to all of the different kinds of supercomputers
in your life, and making plots of y vs x for a given z, but then b vs t for a
given y.

It's like a little robot that does repetitive things for you so you don't get
tenosynovitis.

The Big Ideas of bake:
1 bake can take a file with keys in it and give you a copy with each key
  swapped out for a value
2 bake can take a file, and make different copies of it, with keys replaced
  for different values in each copy
3 bake can do complex things with values by enclosing keys in values
4 bake can do all sorts of operations on each generated set of files
5 bake is easy to extend; things that you would do at the command line can
  easily be rolled into custom bake commands
6 You can make your own version of bake that can do everything that plain bake
  can do, plus things that you can do in Python that would be difficult to do
  through the plain bake command line interface

To learn more about how to use bake, look at the examples in the examples/
directory; in each directory there is an example.txt file that has tutorial
instructions on things to try with that example.
