### Introduction

Bake runs a set of code a number of times, each time having a different set of
input values.

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

***

For a quick example, do:

    $ cd examples/poisson/
    $ bake -m -f bp/sine -e 'python poisson.py; python compare_ideal.py' \
          -o '@label@;sine-n@n@' -o '@n@;3;5;11;21;41'

and you'll see:

    0.233700550136
    0.0530292875455
    0.00826541696623
    0.00205870676453
    0.00051420047815

In one line, bake ran a numerical code for five different test cases, and then
compared those results to the ideal solution, and output those errors in a
list. (Well, the code and the post-processing were already coded up, but bake
coordinated this action five times in one line.)

(If you haven't installed bake yet, instead of `bake`, type 

    $ ../../bake/cmdline.py -m -f bp/sine \
          -e 'python poisson.py; python compare_ideal.py' \
          -o '@label@;sine-n@n@' -o '@n@;3;5;11;21;41'

instead.)

***

### The Big Ideas of bake

1. bake can take a file with keys in it and give you a copy with each key
   swapped out for a value
2. bake can take a file, and make different copies of it, with keys replaced
   for different values in each copy
3. bake can do complex things with values by enclosing keys in values
4. bake can do all sorts of operations on each generated set of files
5. bake is easy to extend; things that you would do at the command line can
   easily be rolled into custom bake commands
6. You can make your own version of bake that can do everything that plain bake
   can do, plus things that you can do in Python that would be difficult to do
   through the plain bake command line interface

To learn more about how to use bake, look at the examples in the examples/
directory; in each directory there is an example.txt file that has tutorial
instructions on things to try with that example. More details on how bake works
are below.

### What you need to run bake

First, you'll need Python.

To install bake, just do:

python setup.py install

You may need to edit your PATH or PYTHONPATH environment variable to find
bake. You can find discussion on this at:

    http://docs.python.org/install/index.html#modifying-python-s-search-path

To set up a project to use bake, you need the
following in the project directory:
* An empty directory named `batch`
* A `bake.cfg` file; you can use the one in the cookie and poisson
  examples. Edit it so that `code_files` indicates the files you'll want bake
  to edit. This list should be comma-separated but with no spaces. It can span
  more than one line, using a `\` at the end of each line.
* A bake parameter (bp) file that describes which parameters should be 
  varied and substituted.

You can then enter keys anywhere in your code files, and, when bake operates,
it will substitute those keys for the values in the bp file.

### How you can help the bake project

Bug reports and code for new features are welcome at this point. However, more
importantly, I have used bake for three projects in my research, but my use of
it has been idiosyncratic. I have tried to separate many of these
idiosyncracies into my own bake front-end that is not included in this project.

What I would appreciate the most are:

1. Questions about why bake could help you
2. Ideas for things you could use bake for
3. Notes on what parts of documentation are unclear

### How bake works

First, bake reads a bp file, line-by-line. If any lines are added through the
overwrite option, those are prepended to the file. Then each line is scanned
(see `bake.mix.parseBPlines`). If the line starts with a `#` (hash) it is a
comment and is ignored. Otherwise, the line should be of the format,

    @key@;value1;value2;...;value_n

The line is split with the semicolons acting as delimiters. The first element,
the key, is put in one list; the remaining items are put in a list that
corresponds to this key. If multiple lines have the same key, only the first of
these lines is read; the rest are ignored. There must be a line with the key
`@label@`.

Bake makes an n-dimensional grid where n is the number of keys with multiple
values. Bake makes an iterator that crawls across this grid, making every
possible combination of values. (See the discussion in
`example/cookie/example.txt` of "Example 4: nuts", in which one key with 3
values and another with 2 lead to 6 jobs.)

For each job, that is, for each location in this grid, a specific set of
values, one for each key, is selected. bake then reads through this list. (See
`bake.mix.TokenValueSubValue()` for this process.) For each value that has a key
in it, the key in that value is substituted for its corresponding value. For
example, if this were a combination of keys and values,

    @spam@: eggs-@foo@-bar
    @foo@: baz

this would reduce to:

    @spam@: eggs-baz-bar
    @foo@: baz

Bake does not care which order these key-value combinations appear in. Also, a
key can have a value that includes a key that has another value including a
key, and so on. Bake iterates over the keys, doing each of these substitutions,
until each value has no keys in it.

This process is aided by having a consistent format for keys. This is specified
by the `pattern` option in `bake.cfg`

By default, bake does these substitutions to generate the list of `@label@`s
for a grid of jobs. When bake does the mix task, it takes these key-value pairs
and substitutes them over all of the source code specified in the `bake.cfg`
`code_files` option.

You can do substitution operations on other strings or even files. An example
of this is in the `bake/examples/poisson/myBake/cmdline.py` file, inside a

    for values in mixIterator:

loop:
        for j in range(0, len(tokens)):
            table_line = table_line.replace(tokens[j], values[j])

This does this key-value substitituion over the string table_line.

#### The `@label@` tag.

For bake to select each job, each job must have a unique label. Job names are
set by the `@label@` tag. This tag should include each key that has multiple
values.

Bake can do operations on each job. The simplest way to do this is with the -e
(execute) task with a shell command as its argument. Bake cd's into each job
directory, and performs the selected command as if you had done this yourself. For example, suppose

    bake -e 'make; ./a.out' -f bp/file

compiles the code and runs it, assuming that you have a makefile and that a.out
is the compiled binary.

Suppose each job outputs a datafile, `foo.txt`. You could then do,

    bake -e 'tail -1 foo.txt' -f bp/file > foolist.txt

and bake would pluck the last line of foo.txt from each job and save this list of last lines of `foo.txt` as `foolist.txt`.

You can assemble your own front-end to bake that can do more sophisticated tasks; an example of this is in `examples/poisson/myBake/cmdline.py`.

### Uses of bake that may be surprising

#### Using overwrite to use one bp file to do many different grids

As a numericist, I frequently have to vary sets of numerical parameters in
concert, say the stiffness of a material and the shear rate of the fluid around
it. I could use a bake file with each value I might want embedded in it, or I
could explore by leaving the bake file unchanged and doing something like

    bake -m -e 'run code' -o '@stiffness@;1;2;3;4' \
        -o '@shear_rate@;100;200;400'

to set up a grid of 12 runs. Or, if I'm adding a new feature and I want to
compare my results with and without that feature, I often do,

    bake -m -e 'run code' -o '@new_feature@;False;True'

#### Multiple keys can depend on a single key; values can be expressions

It's not unusual for me to need to refine many parameters in tandem, for
example, using smaller grids and timesteps. I can do this with something like:

@label@;refinement-study-@fineness@
@fineness@;3;4;5
@dt@;(1/2^@fineness@)
@dx@;(@length@/2^@fineness@)
@length@;42

This generates three jobs, each with half the timestep and grid spacing as the
one before.

Note that `@dt@` and `@dx@` have values that are mathematical expressions. If a
key appears in source code, the code behaves exactly as if it really had the
values substituted in by hand, so as long as these expressions are valid in
whatever language of code they wind up in, that's fine. I recommend enclosing
complicated values in parentheses to preserve order of operations in complex
expressions.

#### Actual restrictions on the value of `@label@`

Earlier, it was stated that each key with multiple values should appear in the
value of the `@label@` key. Technically, this is not exactly true. Each
generated `@label@` should be unique. In this simple bp file:

    @a@;b-@c@
    @c@;d;e

the key, `@c@` has multiple values; this is also captured in `@a@`'s multiple
values, so the label could be:

    @label@;name-@a@

or

    @label@;name-@c@

and bake would be happy. This use is rare, though.
