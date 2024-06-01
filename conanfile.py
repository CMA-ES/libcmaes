import re, os, functools

from conan import ConanFile, tools
from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.files import load

from conan.tools.scm import Git


class CmaesConan(ConanFile):
    name = "libcmaes"

    generators = "CMakeDeps"

    # Optional metadata
    license = "MIT"
    author = "<Put your name here> <And your email here>"
    url = "https://github.com/CMA-ES/libcmaes"
    description = "libcmaes is a multithreaded C++11 library with Python bindings for high performance blackbox stochastic optimization using the CMA-ES algorithm for Covariance Matrix Adaptation Evolution Strategy"
    topics = ("<Put some tag here>", "<here>", "<and here>")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "openmp": [True, False],
        "surrog": [True, False],
        "enable_tests": [True, False],
    }
    default_options = {
        "shared": True,
        "openmp": True,
        "surrog": True,
        "enable_tests": False,
        "boost/*:without_python": False,
    }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "include/*",
        "libcmaes-config.cmake.in",
        "src/*",
        "libcmaes.pc.in",
    )

    def build_requirements(self):
        if self.options.enable_tests:
            self.test_requires("gflags/2.2.2")
            self.test_requires("boost/1.74.0")

    def requirements(self):
        self.requires("eigen/3.4.0", transitive_headers=True)

    def set_version(self):
        content = load(self, os.path.join(self.recipe_folder, "CMakeLists.txt"))
        value = re.search(r"set\(libcmaes_VERSION (.*)\)", content)
        extracted_version = value.group(1).strip()

        is_git_tag = False
        git = Git(self, folder=self.recipe_folder)
        try:
            git.run("describe --exact-match --tags")
            is_git_tag = True
        except Exception:
            is_git_tag = False

        if is_git_tag:
            self.version = extracted_version
        else:
            # if not tag -> pre-release version
            commit_hash = git.get_commit()[:8]
            self.version = f"{extracted_version}.{commit_hash}"

    def config_options(self):
        pass

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["LIBCMAES_BUILD_EXAMPLES"] = False
        tc.variables["LIBCMAES_BUILD_SHARED_LIBS"] = self.options.shared
        tc.variables["LIBCMAES_USE_OPENMP"] = self.options.openmp
        tc.variables["LIBCMAES_ENABLE_SURROG"] = self.options.surrog
        tc.variables["LIBCMAES_BUILD_PYTHON"] = self.options.enable_tests
        tc.variables["LIBCMAES_BUILD_TESTS"] = self.options.enable_tests
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["cmaes"]
        self.cpp_info.set_property("cmake_target_name", "libcmaes::cmaes")
        # self.cpp_info.requires = ["eigen::eigen"]
