This script takes bake parameter files, which are lists of tags and values,
then does something for each combination of values.  Okay, that's pretty vague
and broad.  It's useful, honest, for doing repetitive find-and-replace
operations, wrangling data out of oodles of really-similar-but-sublty-
different-in-two-variables-but-not-the-other-five sets of data, doing the
accounting on submitting jobs to all of the different kinds of supercomputers
in your life, and making plots of y vs x for a given z, but then b vs t for a
given y.

It's like a little robot that does repetitive things for you so you don't get
carpal tunnel.

The Big Ideas of bake:
1 bake can take a file with @tags@ in it and give you a copy with each @tag@
  swapped out for a value
2 bake can take a file, and make different copies of it, with @tags@ replaced
  for different values in each copy
3 bake can do complex things with values by enclosing tags in values
4 bake can do all sorts of operations on each generated set of files
5 bake is easy to extend; things that you would do at the command line can
  easily be rolled into custom bake commands
6 You can make your own version of bake that can do everything that plain bake
  can do, plus things that you can do in Python that would be difficult to do
  at the command line
