#!/usr/bin/env python

from __future__ import print_function
import os, shutil
import os.path as P
import re, sys
from glob import glob
import subprocess
from BinaryBuilder import CMakePackage, GITPackage, Package, stage, warn, \
     PackageError, HelperError, SVNPackage, Apps, write_vw_config, write_asp_config
from BinaryDist import fix_install_paths

def strip_flag(flag, key, env):
    ret = []
    hit = None
    if not key in env:
        return
    for test in env[key].split():
        m = re.search(flag, test)
        if m:
            hit = m
        else:
            ret.append(test)
    if ret:
        env[key] = ' '.join(ret).strip()
    else:
        del env[key]
    return hit, env

class ccache(Package):
    src     = 'http://samba.org/ftp/ccache/ccache-3.1.9.tar.bz2'
    chksum  = 'e80a5cb7301e72f675097246d722505ae56e3cd3'

class m4(Package):
    src     = 'http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.gz'
    chksum  = '4f80aed6d8ae3dacf97a0cb6e989845269e342f0'

class libtool(Package):
    src     = 'http://ftpmirror.gnu.org/libtool/libtool-2.4.2.tar.gz'
    chksum  = '22b71a8b5ce3ad86e1094e7285981cae10e6ff88'

class autoconf(Package):
    src='http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz'
    chksum  = '562471cbcb0dd0fa42a76665acf0dbb68479b78a'

class automake(Package):
    src='ftp://ftp.gnu.org/gnu/automake/automake-1.14.1.tar.gz'
    chksum  = '0bb1714b78d70cab9907d2013082978a28f48a46'

class chrpath(Package):
    src     = 'http://ftp.debian.org/debian/pool/main/c/chrpath/chrpath_0.13.orig.tar.gz'
    chksum  = '11ff3e3dda2acaf1e529475f394f74f2ef7a8204'
    # chrpath pollutes the doc folder
    @stage
    def install(self):
        super(chrpath, self).install()
        cmd = ['rm', '-vrf'] + glob(P.join( self.env['INSTALL_DIR'], 'doc', 'chrpath*' ))
        self.helper(*cmd)

class bzip2(Package):
    src     = 'http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz'
    chksum  = '3f89f861209ce81a6bab1fd1998c0ef311712002'
    def configure(self): pass
    @stage
    def install(self):
        # Copy just the things we need.
        self.helper(*['mkdir','-p',P.join(self.env['INSTALL_DIR'],'include')]);
        self.helper(*['mkdir','-p',P.join(self.env['INSTALL_DIR'],'lib')]);
        self.helper(*['mkdir','-p',P.join(self.env['INSTALL_DIR'],'bin')]);
        cmd = ['cp', '-vf'] + glob(P.join(self.workdir, '*.h')) + \
              [P.join(self.env['INSTALL_DIR'], 'include')]
        self.helper(*cmd)
        cmd = ['cp', '-vf'] + glob(P.join(self.workdir, 'lib*')) + \
              [P.join(self.env['INSTALL_DIR'], 'lib')]
        self.helper(*cmd)
        cmd = ['cp', '-vf', P.join(self.workdir, 'bzip2'),
               P.join(self.env['INSTALL_DIR'], 'bin')]
        self.helper(*cmd)

class pbzip2(Package):
    src     = 'http://compression.ca/pbzip2/pbzip2-1.1.6.tar.gz'
    chksum  = '3b4d0ffa3ac362c3702793cc5d9e61664d468aeb'
    def configure(self): pass
    def compile(self):
        self.helper('sed','-ibak','-e','s# g++# %s#g' % self.env['CXX'],
                    'Makefile');
        cflags = 'CFLAGS = -I' + P.join(self.env['INSTALL_DIR'], 'include') + \
                 ' -L' + P.join(self.env['INSTALL_DIR'], 'lib') + ' '
        self.helper('sed','-ibak','-e','s#CFLAGS = #%s#g' % cflags, 'Makefile')
        super(pbzip2, self).compile()
    def install(self):
        # Copy just the things we need.
        cmd = ['cp', '-vf', P.join(self.workdir, 'pbzip2'),
               P.join(self.env['INSTALL_DIR'], 'bin')]
        self.helper(*cmd)


class tiff(Package):
    src     = 'http://download.osgeo.org/libtiff/tiff-4.0.3.tar.gz'
    chksum  = '652e97b78f1444237a82cbcfe014310e776eb6f0'

    def configure(self):
        super(tiff, self).configure(
            with_ = ['jpeg', 'png', 'zlib'],
            without = ['x'],
            enable=('shared',),
            disable = ['static', 'lzma', 'cxx', 'logluv'])


class qt(Package):
    src     = 'http://download.qt-project.org/official_releases/qt/4.8/4.8.6/qt-everywhere-opensource-src-4.8.6.tar.gz'
    chksum  = 'ddf9c20ca8309a116e0466c42984238009525da6' #SHA-1 Hash
    patches = 'patches/qt'
    patch_level = '-p0'

    def __init__(self, env):
        super(qt, self).__init__(env)

        # Qt can only be built on OSX with an Apple Compiler. If the
        # user overwrote the compiler choice, we must revert here. The
        # problem is -fconstant-cfstrings. Macports also gives up in
        # this situation and blacks lists all Macport built compilers.
        if self.arch.os == 'osx':
            self.env['CXX']='c++'
            self.env['CC']='cc'

    @stage
    def configure(self):
        # The default confs override our compiler choices.
        self.helper('sed','-ibak','-e','s# g++# %s#g' % self.env['CXX'], '-e', 's# gcc# %s#g' % self.env['CC'], 'mkspecs/common/g++-base.conf')
        cmd = './configure -opensource -fast -confirm-license -nomake demos -nomake examples -nomake docs -nomake translations -no-webkit -prefix %(INSTALL_DIR)s -no-script -no-scripttools -no-openssl -no-libjpeg -no-libmng -no-libpng -no-libtiff -no-cups -no-nis -no-opengl -no-openvg -no-phonon -no-phonon-backend -no-sql-psql -no-dbus' % self.env
        args = cmd.split()
        if self.arch.os == 'osx':
            args.append('-no-framework')
            args.extend(['-arch',self.env['OSX_ARCH']])
        self.helper(*args)

    @stage
    def install(self):
        # Call the install itself afterward
        super(qt, self).install()

class zlib(Package):
    src     = 'http://downloads.sourceforge.net/libpng/zlib-1.2.8.tar.gz'
    chksum  = 'a4d316c404ff54ca545ea71a27af7dbc29817088'

    @stage
    def configure(self):
        super(zlib,self).configure(other=('--shared',))

    @stage
    def install(self):
        super(zlib, self).install()
        self.helper(*['rm', P.join(self.env['INSTALL_DIR'], 'lib', 'libz.a')])

class jpeg(Package):
    src     = 'http://www.ijg.org/files/jpegsrc.v8d.tar.gz'
    chksum  = 'f080b2fffc7581f7d19b968092ba9ebc234556ff'
    patches = 'patches/jpeg8'

    def configure(self):
        super(jpeg, self).configure(enable=('shared',), disable=('static',))

class png(Package):
    src    = 'http://downloads.sourceforge.net/libpng/libpng-1.6.7.tar.gz'
    chksum = '22fcd1aaab3d8f4b98f43e5b301cc4fd7cc15722'

    def configure(self):
        super(png,self).configure(disable='static')


class polyview(GITPackage):
    src     = 'https://github.com/oleg-alexandrov/PolyView.git'

    def configure(self):
        cmd = ('qmake', 'polyview.pro')
        self.helper(*cmd)

    @stage
    def compile(self, cwd=None):
        super(polyview, self).compile(cwd)
