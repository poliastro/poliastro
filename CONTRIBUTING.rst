Contributing
============

poliastro is a community project, hence all contributions are more than
welcome!

Bug reporting
-------------

Not only things break all the time, but also different people have different
use cases for the project. If you find anything that doesn't work as expected
or have suggestions, please refer to the `issue tracker`_ on GitHub.

.. _`issue tracker`: https://github.com/poliastro/poliastro/issues

Documentation
-------------

Documentation can always be improved and made easier to understand for
newcomers. The docs are stored in text files under the `docs/source`
directory, so if you think anything can be improved there please edit the
files and proceed in the same way as with `code writing`_.

The Python classes and methods also feature inline docs: if you detect
any inconsistency or opportunity for improvement, you can edit those too.

Besides, the `wiki`_ is open for everybody to edit, so feel free to add
new content.

Code writing
------------

.. image:: https://img.shields.io/waffle/label/poliastro/poliastro/1%20-%20Ready.svg?style=flat-square
   :target: https://waffle.io/poliastro/poliastro
   :alt: 'Stories in Ready'

Code contributions are welcome! If you are looking for a place to start,
help us fixing bugs in poliastro and check out the `"easy" tag`_. Those
should be easier to fix than the others and require less knowledge about the
library.

.. _`"easy" tag`: https://github.com/poliastro/poliastro/issues?q=is%3Aissue+is%3Aopen+label%3Aeasy

If you are hesitant on what IDE or editor to use, just choose one that
you find comfortable and stick to it while you are learning. People have
strong opinions on which editor is better so I recommend you to ignore
the crowd for the time being - again, choose one that you like :)

If you ask me for a recommendation, I would suggest PyCharm (complete
IDE, free and gratis, RAM-hungry) or vim (powerful editor, very lightweight,
steep learning curve). Other people use Spyder, emacs, gedit, Notepad++,
Sublime, Atom...

You will also need to understand how git works. git is a decentralized
version control system that preserves the history of the software, helps
tracking changes and allows for multiple versions of the code to exist
at the same time. If you are new to git and version control, I recommend
following `the Try Git tutorial`_.

.. _`the Try Git tutorial`: https://try.github.io/

If you already know how all this works and would like to contribute new
features then that's awesome! Before rushing out though please make sure it
is within the scope of the library so you don't waste your time -
`the mailing list`_ or `the chat`_ are good places to ask.

.. _`the mailing list`: https://groups.io/g/poliastro-dev
.. _`the chat`: https://riot.im/app/#/room/#poliastro:matrix.org

If the feature you suggest happens to be useful and within scope, you will
probably be advised to create a new `wiki`_ page with some information
about what problem you are trying to solve, how do you plan to do it and
a sketch or idea of how the API is going to look like. You can go there
to read other good examples on how to do it. The purpose is to describe
without too much code what you are trying to accomplish to justify the
effort and to make it understandable to a broader audience.

.. _`wiki`: https://github.com/poliastro/poliastro/wiki

All new features should be thoroughly tested, and in the ideal case the
coverage rate should increase or stay the same. Automatic services will ensure
your code works on all the operative systems and package combinations
poliastro support - specifically, note that poliastro is a Python 3 only
library.

Development environment
-----------------------

These are some succint steps to set up a development environment:

1. `Install git <https://git-scm.com/>`_ on your computer.
2. `Register to GitHub <https://github.com/>`_.
3. `Fork poliastro <https://help.github.com/articles/fork-a-repo/>`_.
4. `Clone your fork <https://help.github.com/articles/cloning-a-repository/>`_.
5. Install it in development mode using
   :code:`pip install --editable /path/to/poliastro/[dev]` (this means that the
   installed code will change as soon as you change it in the download
   location).
6. Create a new branch.
7. Make changes and commit.
8. `Push to your fork <https://help.github.com/articles/pushing-to-a-remote/>`_.
9. `Open a pull request! <https://help.github.com/articles/creating-a-pull-request/>`_

For more detailed explanations, please check out the `Astropy development docs`__.

.. __: http://docs.astropy.org/en/stable/development/workflow/development_workflow.html
