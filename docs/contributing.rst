============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
^^^^^^^^^^^

Report bugs at https://github.com/GeoscienceAustralia/PyRate/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

For an example of how to report a bug please see: https://github.com/GeoscienceAustralia/PyRate/issues/146

Fix Bugs
^^^^^^^^

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
^^^^^^^^^^^^^^^^^^

Look through the GitHub issues for new features. Anything tagged with
"enhancement" is open to whoever wants to implement it.

Write Documentation
^^^^^^^^^^^^^^^^^^^

PyRate could always use more documentation, whether as part of the
official PyRate docs, in docstrings, or even on the web in blog posts,
articles etc.

Submit Feedback
^^^^^^^^^^^^^^^

The best way to send feedback is to file an Issue_.

.. _Issue: https://github.com/GeoscienceAustralia/PyRate/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `PyRate` for local development.

1. Fork_ the `PyRate` repository on GitHub.
2. Clone your fork locally

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git

3. Create a branch for local development

::

    git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

4. When you're finished making changes, check that your changes pass style and unit
   tests. A suite of tests have been developed for use in testing PyRate functionality
   and for further code development. The tests use `pytest <http://doc.pytest.org/en/latest/>`__
   and can be found in the *tests/* directory. A small test dataset is included in the *tests/test\_data/* directory.

::

    cd PyRate
    pytest tests/
    pip install tox
    tox

5. Commit your changes and push your branch to GitHub

::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

.. _Fork: https://help.github.com/articles/fork-a-repo/

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the documentation should be updated.
   Put your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.3+ and for PyPy.

   Check https://travis-ci.org/GeoscienceAustralia/PyRate
   under pull requests for active pull requests or run the ``tox`` command and
   make sure that the tests pass for all supported Python versions.
