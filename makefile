###############################################################################
#                                                                             #
#  Copyright © 2017-2018 Gustavo Thebit Pfeiffer / Y. Sato laboratory         #
#                                                                             #
#                                                                             #
#  This file is part of PS-MCLS2018.                                          #
#                                                                             #
#  PS-MCLS2018 is free software: you can redistribute it and/or modify it     #
#  under the terms of the GNU Lesser General Public License as published by   #
#  the Free Software Foundation, either version 3 of the License, or (at your #
#  option) any later version.                                                 #
#                                                                             #
#  PS-MCLS2018 is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY #
#  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public     #
#  License for more details.                                                  #
#                                                                             #
#  You should have received a copy of the GNU Lesser General Public License   #
#  along with PS-MCLS2018. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################

v = main.cpp
f = out-$(v)
all:
	-rm -r $(f)
	mkdir $(f)
	cp $(v) core.h methods.h problems.h $(f)
	g++ $(v) plotter/plotter.cpp -O3 -std=c++11 -lSDL -lSDL_ttf -o main
	(cd $(f); ../main;)
	-inkscape $(f)/out.svg -E $(f)/out.eps --export-ignore-filters --export-ps-level=3
	-epstopdf $(f)/out.eps
	-rm $(f)/out.eps

install:
	apt-get install g++ libsdl1.2-dev libsdl-ttf2.0-dev libeigen3-dev inkscape texlive-font-utils


