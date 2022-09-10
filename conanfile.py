import re, os, functools

from conans import tools as tools
from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conans.tools import load

class CmaesConan(ConanFile):
    name = "libcmaes"

    generators = "CMakeDeps", "CMakeToolchain"

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
        "surrog": [True, False]
        }
    default_options = {
        "shared": True, 
        "openmp": True,
        "surrog": True
        }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "cmake/*", "include/*", "libcmaes-config.cmake.in", "src/*", "libcmaes.pc.in"

    def build_requirements(self):
        pass

    def requirements(self):
        self.requires("eigen/3.4.0")
        

    def set_version(self):
        content = load(os.path.join(self.recipe_folder, "CMakeLists.txt"))
        value=re.search(r"set\(libcmaes_VERSION (.*)\)", content)
        print(value)
        extracted_version  = value.group(1).strip()
        
        git = tools.Git(folder=self.recipe_folder)
        if (git.get_tag() != None):
            # depending on your workflow you could also
            # set a non-beta version when on main branch or
            # on a certain release branch
            self.version = extracted_version
        else:
            # if not tag -> pre-release version
            commit_hash = git.get_commit()[:8]
            branch_name = git.get_branch()[:9]
            self.version = f"{extracted_version}-{branch_name}.{commit_hash}"

    def config_options(self):
        pass

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables['LIBCMAES_BUILD_EXAMPLES']=False
        tc.variables['LIBCMAES_BUILD_SHARED_LIBS']= self.options.shared
        tc.variables['LIBCMAES_USE_OPENMP'] = self.options.openmp
        tc.variables['LIBCMAES_ENABLE_SURROG'] = self.options.surrog
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
        self.cpp_info.set_property("cmake_target_aliases",["libcmaes::cmaes"])
        self.cpp_info.requires = ["eigen::eigen"]
