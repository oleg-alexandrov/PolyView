#!/bin/sh

cp -fv gui/polyview debian/usr/bin/
dpkg-deb --build debian
mv -fv debian.deb polyview_0.5.deb

# Notes to self:

# Zip the debian directory:
# tar czfv debian.tgz debian

# Install (need to be root)
# dpkg -i polyView_0.5.deb

# Uninstall
# dpkg -r polyview

# See what libraries an executable depends on:
# ldd exeName

# See which libraries are installed in Ubuntu:
# dpkg --get-selections

# Find the Ubuntu version:
# cat /etc/issue


