# Contributing

poliastro is a community project, all contributions are more than welcome!

## Bug reporting

Not only things can break, but also different people have
different use cases for the project. If you find anything that doesn't
work as expected or have suggestions, open a new issue on our
[issue tracker](https://github.com/poliastro/poliastro/issues).

## Documentation

Documentation can always be expanded and improved. The docs are stored in
text files under the `docs/source` directory, so if you have any suggestions,
edit the files and proceed in the same way as with [code writing](#code-writing).

The Python classes and methods also feature inline docs in the form of
docstrings. If you detect any inconsistencies or opportunities for
improvement, you can edit those too.

To build the docs, you must first create a development environment (see
below) and then in the `docs/` directory run:

```console
$ cd docs
$ make html
```

After this, the new docs will be inside `build/html`. You can open them
by running an HTTP server:

```console
$ cd build/html
$ python -m http.server
Serving HTTP on 0.0.0.0 port 8000 ...
```

And point your browser to <http://0.0.0.0:8000>.

## Small scripts

We would love to give your Astrodynamics scripts a home!
Please head to [our `contrib/` directory](https://github.com/poliastro/poliastro/tree/main/contrib)
for further information.

## Code writing

Code contributions are welcome! If you are looking for a place to start,
check out the ["good-first-issue" label](https://github.com/poliastro/poliastro/labels/good%20first%20issue)
on our issue tracker. Those tasks should be easier to fix than the others
and require less knowledge about the library.

If you are hesitant on what IDE or editor to use, choose one that
you find comfortable and stick to it while you are learning.
If you want a recommendation, have a look at PyCharm, Visual Studio Code, or Spyder.

You will also need to understand how git works. git is a decentralized
version control system that preserves the history of the software, helps
tracking changes and allows for multiple versions of the code to exist
at the same time. If you are new to git and version control, have a look at
[this Getting Started guide](https://docs.github.com/en/get-started/getting-started-with-git).

If you already know how all this works and would like to contribute new
features then that's awesome! Before rushing out though make
sure it is within the scope of the library so you don't waste your time
 - [the mailing list](https://groups.io/g/poliastro-dev) or [the
chat](http://chat.poliastro.space/) are good places to ask.

All new features should be thoroughly tested, and in the ideal case the
coverage rate should increase or stay the same. Automatic services will
ensure your code works on all the operative systems and package
combinations poliastro support.

## Development environment

These are some succint steps to set up a development environment:

1. [Install git](https://git-scm.com/) on your computer.
2. [Sign up to GitHub](https://github.com/).
3. [Fork poliastro](https://help.github.com/articles/fork-a-repo/).
4. [Clone your fork](https://help.github.com/articles/cloning-a-repository/).
5. Create a Python virtual environment using `python -m venv .venv`
6. Activate it using `source .venv/bin/activate`
7. Upgrade the development dependencies using `python pip install -U pip setuptools wheel flit tox`
8. Install the code in development mode using `flit install --simlink`
   (this means that the installed code will change as soon as you change it in the
   download location).
9. Create a new branch using `git switch --create descriptive-name-of-my-change`.
10. Make your changes!
11. Run `tox -e reformat` to make the format consistent.
12. Run `tox -e check` to check all the formatting is right.
13. Commit your changes using `git commit -m 'Add new awesome feature'`.
14. [Push to your fork](https://help.github.com/articles/pushing-to-a-remote/).
15. [Open a pull request!](https://help.github.com/articles/creating-a-pull-request/)

For more detailed explanations, check out the [Astropy development
docs](http://docs.astropy.org/en/stable/development/workflow/development_workflow.html).
