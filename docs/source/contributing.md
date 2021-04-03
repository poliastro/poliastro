# Contribute

poliastro is a community project, hence all contributions are more than
welcome!

## Bug reporting

Not only things break all the time, but also different people have
different use cases for the project. If you find anything that doesn\'t
work as expected or have suggestions, please refer to the [issue
tracker](https://github.com/poliastro/poliastro/issues) on GitHub.

## Documentation

Documentation can always be improved and made easier to understand for
newcomers. The docs are stored in text files under the
*docs/source* directory, so if you think anything can be
improved there please edit the files and proceed in the same way as with
[code writing](#code-writing).

The Python classes and methods also feature inline docs: if you detect
any inconsistency or opportunity for improvement, you can edit those
too.

Besides, the [wiki](https://github.com/poliastro/poliastro/wiki) is open
for everybody to edit, so feel free to add new content.

To build the docs, you must first create a development environment (see
below) and then in the `docs/` directory run:

```bash
$ cd docs
$ make html
```

After this, the new docs will be inside `build/html`. You can open them
by running an HTTP server:

```bash
$ cd build/html
$ python -m http.server 8000 --bind 127.0.0.1
Serving HTTP on 127.0.0.1 port 8000 (http://127.0.0.1:8000/) ...
```

And point your browser to <http://127.0.0.1:8000>.

## Code writing

Code contributions are welcome! If you are looking for a place to start,
help us fixing bugs in poliastro and check out the [\"good first issue\"
tag](https://github.com/poliastro/poliastro/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22).
Those should be easier to fix than the others and require less knowledge
about the library.

If you are hesitant on what IDE or editor to use, just choose one that
you find comfortable and stick to it while you are learning. People have
strong opinions on which editor is better so I recommend you to ignore
the crowd for the time being - again, choose one that you like :)

If you ask me for a recommendation, I would suggest PyCharm (complete
IDE, free and gratis, RAM-hungry) or vim (powerful editor, very
lightweight, steep learning curve). Other people use Spyder, emacs,
gedit, Notepad++, Sublime, Atom\...

You will also need to understand how git works. git is a decentralized
version control system that preserves the history of the software, helps
tracking changes and allows for multiple versions of the code to exist
at the same time. If you are new to git and version control, I recommend
following [the Try Git tutorial](https://try.github.io/) and [Git
Workflow tutorial](https://www.asmeurer.com/git-workflow/).

If you already know how all this works and would like to contribute new
features then that\'s awesome! Before rushing out though please make
sure it is within the scope of the library so you don\'t waste your time
-[the mailing list](https://groups.io/g/poliastro-dev) or [the
chat](http://chat.poliastro.space/) are good places to ask.

If the feature you suggest happens to be useful and within scope, you
will probably be advised to create a new
[wiki](https://github.com/poliastro/poliastro/wiki) page with some
information about what problem you are trying to solve, how do you plan
to do it and a sketch or idea of how the API is going to look like. You
can go there to read other good examples on how to do it. The purpose is
to describe without too much code what you are trying to accomplish to
justify the effort and to make it understandable to a broader audience.

All new features should be thoroughly tested, and in the ideal case the
coverage rate should increase or stay the same. Automatic services will
ensure your code works on all the operative systems and package
combinations poliastro support - specifically, note that poliastro is a
Python 3 only library.

## Development environment

These are some succinct steps to set up a development environment:

1.  [Install git](https://git-scm.com/) on your computer.
2.  [Register to GitHub](https://github.com/).
3.  [Fork poliastro](https://help.github.com/articles/fork-a-repo/).
4.  [Clone your
    fork](https://help.github.com/articles/cloning-a-repository/).
5.  Install [flit](https://github.com/takluyver/flit) using
    `python -m pip install flit wheel`
6.  Change into the local checkout of your repository and install it in
    development mode using `flit install --symlink` (`--symlink` means
    that the installed code will change as soon as you change it in the
    download location).
7.  Create a new branch.
8.  Make changes and commit.
9.  [Push to your
    fork](https://help.github.com/articles/pushing-to-a-remote/).
10. [Open a pull
    request!](https://help.github.com/articles/creating-a-pull-request/)

For more detailed explanations, please check out the [Astropy
development
docs](https://docs.astropy.org/en/stable/development/workflow/development_workflow.html).
