---@diagnostic disable: undefined-global, undefined-field

add_repositories("my-repo https://github.com/45degree/personal-xmake-repo.git main")

add_rules("mode.debug", "mode.release")
add_requires("openmesh", {debug = true})
add_requires("libigl", "tbb")
add_requires("gtest")
add_requires("progressbar v2.1")

target("poly2tri")
    set_kind("static")
    add_files("poly2tri/poly2tri/**.cc")
    add_includedirs("poly2tri/poly2tri/", { public = true })
target_end()

target("Maps")
    set_kind("static")
    set_languages("c99", "c++17")
    add_files("src/*.cc")
    add_packages("openmesh", "libigl", "tbb", "progressbar")
    add_includedirs("include/", { public = true })
    if is_plat("linux") then
        add_syslinks("pthread")
    end
    add_deps("poly2tri")

    if is_mode("debug") and is_plat("linux") then
        add_cxflags("-pg")
    end

    on_load(function(target)
        os.cp("$(projectdir)/model", target:targetdir())
    end)
target_end()

target("MapsTest")
    set_kind("binary")
    set_languages("c99", "c++17")
    add_files("test/*.cc")
    add_deps("Maps")
    add_packages("gtest", "openmesh", "libigl")
target_end()

target("MapsExample")
    set_kind("phony")
    add_deps("MapsExample1", "MapsExample2")
target_end()

target("MapsExample2")
    set_kind("binary")
    set_languages("c99", "c++17")
    add_files("example/example2.cc")
    add_deps("Maps")
    add_packages("openmesh")
target_end()

target("MapsExample1")
    set_kind("binary")
    set_languages("c99", "c++17")
    add_files("example/example1.cc")
    add_deps("Maps")
    add_packages("openmesh")
target_end()

--
-- If you want to known more usage about xmake, please see https://xmake.io
--
-- ## FAQ
--
-- You can enter the project directory firstly before building project.
--
--   $ cd projectdir
--
-- 1. How to build project?
--
--   $ xmake
--
-- 2. How to configure project?
--
--   $ xmake f -p [macosx|linux|iphoneos ..] -a [x86_64|i386|arm64 ..] -m [debug|release]
--
-- 3. Where is the build output directory?
--
--   The default output directory is `./build` and you can configure the output directory.
--
--   $ xmake f -o outputdir
--   $ xmake
--
-- 4. How to run and debug target after building project?
--
--   $ xmake run [targetname]
--   $ xmake run -d [targetname]
--
-- 5. How to install target to the system directory or other output directory?
--
--   $ xmake install
--   $ xmake install -o installdir
--
-- 6. Add some frequently-used compilation flags in xmake.lua
--
-- @code
--    -- add debug and release modes
--    add_rules("mode.debug", "mode.release")
--
--    -- add macro defination
--    add_defines("NDEBUG", "_GNU_SOURCE=1")
--
--    -- set warning all as error
--    set_warnings("all", "error")
--
--    -- set language: c99, c++11
--    set_languages("c99", "c++11")
--
--    -- set optimization: none, faster, fastest, smallest
--    set_optimize("fastest")
--
--    -- add include search directories
--    add_includedirs("/usr/include", "/usr/local/include")
--
--    -- add link libraries and search directories
--    add_links("tbox")
--    add_linkdirs("/usr/local/lib", "/usr/lib")
--
--    -- add system link libraries
--    add_syslinks("z", "pthread")
--
--    -- add compilation and link flags
--    add_cxflags("-stdnolib", "-fno-strict-aliasing")
--    add_ldflags("-L/usr/local/lib", "-lpthread", {force = true})
--
-- @endcode
--

