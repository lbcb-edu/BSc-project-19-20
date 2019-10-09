# BSc project (computer science - 2018/2019)

Software Design Project is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the fifth semester of the undergraduate study. The main focus is to promote cooperation between students while solving specific problems. Under the supervision of prof. Mile Šikić, students will get familiar with C++, basics of compilation methods, version control, unit tests and continuous integration, and will be introduced to algorithms used in bioinformatics. This project will be executed in several steps, each with defined outcomes which are required for succeeding steps. Instructions and guidelines for the project will be written in this README file which will be updated through the semester.

## Preliminaries

Students are required to get through the following tutorials: [C++](http://www.cplusplus.com/doc/tutorial/), [GitHub](http://rogerdudler.github.io/git-guide/), [CMake](https://cmake.org/cmake-tutorial/), [Googletest](https://github.com/google/googletest/blob/master/googletest/docs/primer.md) and [TravisCI](https://docs.travis-ci.com/user/getting-started/). While writing C++ code, it is advised to follow the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html).

Students will be assigned to one of five teams which are **blue**, **brown**, **orange** and **pink**. Each team will have a separate branch and only team members will have write permission. Although, additional branches can be created if needed, but should have names starting with the team name (e.g. `blue_feature_one`).

## Objective

At the end of the project, students will have implemented several libraries which will enable alignment of a large amount of substrings (of various sizes) to a much larger string from which they originate. The objective is to join these libraries into a single program, often called mapper, in order to map long erroneous fragments from third generation of sequencing technologies to a reference genome, which has various use cases in bioinformatics. A visual example can be seen bellow.

![](misc/sample_mappings.png)
