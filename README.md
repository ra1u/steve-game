
# About

Solution to [blog post](https://rahne.si/optimisation/2026/01/07/steve-ballmer-interview.html) for adversarial binary search game

# How to use

~~~bash
uv run python game.py -h

usage: Steve game [-h] [-n COUNT] [-f FILE_EXPORT] [-b] [-c] [-s]

Solves adverserial bisection search game

options:
  -h, --help            show this help message and exit
  -n, --count COUNT     provide problem size, default = 3
  -f, --file-export FILE_EXPORT
                        export problem in provided mps file - combine with flags -disable-steve/--disable-candidate
  -b, --disable-steve   when enabled, model for Steve strategy is not generated
  -c, --disable-candidate
                        when enabled, model for candidate strategy is not generated
  -s, --dont-solve      when enabled, solver does not run, can be used to just export model in mps file
~~~


Naive way for size 100 (should take approx 1h on laptop)

~~~bash
uv run python game.py -n 100
~~~

To solve as rational problem 

Export problem as mps file

~~~bash
uv run python game.py --file-export "problem.mps" --disable-candidate --dont-solve --count 100
~~~

solve with [SoPlex](https://github.com/scipopt/soplex/tree/master)

~~~bash
soplex  --loadset=exact.set  -X problem.mps
~~~

Requires setting file [exact.set](https://github.com/scipopt/soplex/blob/master/settings/exact.set)  for forcing into rational solver mode.