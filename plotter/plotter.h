//Copyright © 2017-2018 Gustavo Thebit Pfeiffer / Y. Sato laboratory 
/*
    This file is part of PS-MCLS2018.

    PS-MCLS2018 is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    PS-MCLS2018 is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PS-MCLS2018. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>
#include <string>
#include <sstream>
#include <list>

class Plotter {
public:
	int margin_left; //default: 90
	int margin_right; //default: 300
	int margin_top; //default: 30
	int margin_bottom; //default: 60

	double xmin0, ymin0; //+infinite by default
	double xmax0, ymax0; //-infinite by default
	
	bool quit; //Plotter::refresh will set 'quit' to true if ESC was pressed or the window was closed.
	
	std::string xaxisname; //x axis name
	std::string yaxisname; //y axis name
	std::string title; //graph title

	struct Curve {
		double *x; //x values. When left NULL, incremental values are used (1, 2, 3, ...)
		double *y; //y values. Always allocate.
		double *y_u; //upper error (use y_u > y). Not drawn when left NULL.
		double *y_d; //lower error (use y_d < y). Not drawn when left NULL.
		int n; //size of the vectors above (more precisely, this is not the actual allocated memory size, but the number of points available for rendering)
		Uint32 col; //color of the main line. Black by default.
		Uint32 col_u; //color of the upper error line. Grey by default.
		Uint32 col_d; //color of the lower error line. Grey by default.
		
		std::string name; //curve name for legend
		
		Curve() {
			x = NULL;
			y = NULL;
			y_u = NULL;
			y_d = NULL;
			n = 0;
			col = 0;
			col_u = col_d = 0x808080;
		};
	};

private:

	SDL_Surface * srf;
	TTF_Font * font;
	::std::ostringstream tosvg;
	bool tosvg_on;

	std::list<Curve> curves;
	
	double x_a;
	double x_b;
	double y_a;
	double y_b;

	void safe_draw_line(SDL_Surface* srf, double xx1, double yy1, double xx2, double yy2, Uint32 col = 0);
	void draw_line_crop(SDL_Surface * srf, double x1, double y1, double x2, double y2, Uint32 col, double limleft, double limright, double limup, double limdown);

	void drawText(int x, int y, const char * text, int mode = 0, int fontsize = 12);

	void drawVer(double x);
	void drawHor(double y);
	
	void draw_main(const Curve & c);
	void draw_ud(const Curve & c);

	void update_ab();
	void clearGraph();
	void drawEverything();

public:
	//Call refresh to draw (update) graph and handle user interface events
	//If not called, graph will never be updated
	void refresh();

	//Adds a new curve to the list of curves and returns its pointer.
	Curve * addCurve();

	//Saves graph as SVG
	//exportSVG is called automatically when the window is closed
	void exportSVG();
	std::string svg_filename; //"out.svg" by default
	bool svg_high_quality; //when true, SVG graph coordinates are exported in double precision. False by default.

	//Input: graph width/height in pixels
	Plotter(int width = 1000, int height = 600);

	virtual ~Plotter();

};

