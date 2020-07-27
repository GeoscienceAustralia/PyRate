============
Contributing
============

Contributions to `PyRate` are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Types of Contributions
----------------------

You can contribute to `PyRate` in many ways:

Report Bugs
^^^^^^^^^^^

Report bugs on the `Github issues`_ page. If you are reporting a bug, please include:

.. _`Github issues`: https://github.com/GeoscienceAustralia/PyRate/issues

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

For an example of how to report a bug please see: https://github.com/GeoscienceAustralia/PyRate/issues/146

Fix Bugs
^^^^^^^^

Look for any `GitHub issues`_ with the "bug" label.

Implement Features
^^^^^^^^^^^^^^^^^^

Look for any `GitHub issues`_ with the "enhancement" label.

Write Documentation
^^^^^^^^^^^^^^^^^^^

`PyRate` could always use more documentation, whether as part of the
official documentation, in docstrings, or even on the web in blog posts,
articles etc.

Submit Feedback
^^^^^^^^^^^^^^^

General feedback or questions about `PyRate` can be emailed to the
`Geoscience Australia InSAR Team`_.

.. _`Geoscience Australia InSAR Team`: mailto:insar@ga.gov.au

The best way to send feedback about a bug in the software is to lodge an Issue_.

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
2. Clone your fork locally:

::

    cd ~
    git clone git@github.com:<account_where_you_forked>/PyRate.git

3. Create a branch for local development

::

    git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

4. When you have finished making changes, check that your changes pass against unit
   tests. A suite of tests have been developed for use in testing `PyRate` functionality
   and for further code development. The tests use `pytest <http://doc.pytest.org/en/latest/>`__
   and can be found in the ``PyRate/tests/`` directory. A small test dataset used by
   the test suite is included in the repository in the ``PyRate/tests/test_data/small_test`` directory.

::

    cd ~/PyRate
    # 'not slow' avoids running those tests marked as 'slow'
    pytest tests/ -m "not slow"

5. If the tests pass, commit your changes and push your branch to GitHub:

::

    git add .
    git commit -m "Short description of your changes."
    git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

.. _Fork: https://help.github.com/articles/fork-a-repo/

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests that cover new functionality.
2. If the pull request adds functionality, the documentation should be updated.
   Put your new functionality into a function with a docstring.
3. The pull request should work for Python 3.6+.

   Check https://travis-ci.org/GeoscienceAustralia/PyRate
   under pull requests for active pull requests or run the ``tox`` command and
   make sure that the tests pass for all supported Python versions.
