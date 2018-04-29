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

#include "plotter.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
using namespace std;

#define FONTNAME "DejaVuSans"
#define FONTPATH "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

const char * HEXA[] = {"0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"};

template <typename T>
std::string to_stdstring(const T & x) {
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

template <typename T>
inline T sqr(const T & x) {
	return x*x;
}

void draw_line(SDL_Surface* srf,double xx1, double yy1, double xx2, double yy2, Uint32 linecol = 0)
{
	double xl = xx2-xx1;
	double yl = yy2-yy1;
	double imax = sqrt(sqr(xl)+sqr(yl));
	xl /= imax;
	yl /= imax;
	double x = xx1;
	double y = yy1;
	for(int i = 0; i < imax; i++) {
		((Uint32*)(srf->pixels))[srf->w*(int)y+(int)x] = linecol;
		x += xl;
		y += yl;
	}
}

inline std::string to_hexa(const Uint32 & col) {
	return std::string("") +
		HEXA[(col/0x100000)%16] + 
		HEXA[(col/0x10000)%16] + 
		HEXA[(col/0x1000)%16] + 
		HEXA[(col/0x100)%16] + 
		HEXA[(col/0x10)%16] + 
		HEXA[col%16];
}
 
void Plotter::safe_draw_line(SDL_Surface* srf, double xx1, double yy1, double xx2, double yy2, Uint32 col)
{
    if(xx1 >= 0 &&
        yy1 >= 0 &&
        xx1 < srf->w &&
        yy1 < srf->h &&
        xx2 >= 0 &&
        yy2 >= 0 &&
        xx2 < srf->w &&
        yy2 < srf->h) {
        if(tosvg_on) {
			tosvg << 
				"  <line x1=\"" + 
				to_stdstring(xx1) + 
				"\" y1=\"" + 
				to_stdstring(yy1) + 
				"\" x2=\"" + 
				to_stdstring(xx2) + 
				"\" y2=\"" + 
				to_stdstring(yy2) + 
				"\" stroke=\"#" + 
				to_hexa(col) + 
				"\" stroke-width=\"1\" />\n";
        }
	    draw_line(srf, (int)xx1, (int)yy1, (int)xx2, (int)yy2, col);
    }
}

inline void Plotter::draw_line_crop(SDL_Surface * srf, double x1, double y1, double x2, double y2, Uint32 col, double limleft, double limright, double limup, double limdown) {
	if (x1 >= limleft && x2 >= limleft &&
	    x1 <= limright && x2 <= limright &&
	    y1 >= limup && y2 >= limup &&
	    y1 <= limdown && y2 <= limdown) {
		if(x1 >= 0 &&
		    y1 >= 0 &&
		    x1 < srf->w &&
		    y1 < srf->h &&
		    x2 >= 0 &&
		    y2 >= 0 &&
		    x2 < srf->w &&
		    y2 < srf->h) {
			draw_line(srf, (int)x1, (int)y1, (int)x2, (int)y2, col);
		}
    }
}

void Plotter::draw_ud(const Curve & c) {

	double limright = srf->w - margin_right;
	double limdown = srf->h - margin_bottom;

	if(c.y_u && c.y_d) {
		{
			double xl = x_a * (c.x?c.x[0]:1) + x_b;
			double yul = y_a * c.y_u[0] + y_b;
			double ydl = y_a * c.y_d[0] + y_b;
			for(int i = 1; i < c.n; i++) {
				double x = x_a * (c.x?c.x[i]:(i+1)) + x_b;
				double yu = y_a * c.y_u[i] + y_b;
				double yd = y_a * c.y_d[i] + y_b;
				draw_line_crop(srf, xl, yul, x, yu, c.col_u, margin_left, limright, margin_top, limdown);
				draw_line_crop(srf, xl, ydl, x, yd, c.col_d, margin_left, limright, margin_top, limdown);
				xl = x;
				yul = yu;
				ydl = yd;
			}
		}

		if(tosvg_on)
		{
			tosvg << "<polygon points=\"";
			double limup = margin_top;
			double limleft = margin_left;
			
			double xl = 1.0/0.0;

			int imax = c.n;
			for(int i = 0; i < c.n; i++) {
				double x = x_a * (c.x?c.x[i]:(i+1)) + x_b;
				double y = y_a * c.y_u[i] + y_b;
				if(!(y <= limdown))
					y = limdown;
				if(!(y >= limup))
					y = limup;
				
				if(!(c.y[i]-c.y[i]==0)) { //check for NaN or inf
					imax = i;
					break;
				}

				if(x >= 0 &&
					y >= 0 &&
					x < srf->w &&
					y < srf->h &&
					!((int)x == (int)xl)) {
					tosvg << x << "," << y << " ";
				}
				xl = x;
			}

			xl = 1.0/0.0;

			for(int i = imax-1; i >= 0; i--) {
				double x = x_a * (c.x?c.x[i]:(i+1)) + x_b;
				double y = y_a * c.y_d[i] + y_b;
				if(!(y <= limdown))
					y = limdown;
				if(!(y >= limup))
					y = limup;
				
				if(x >= 0 &&
					y >= 0 &&
					x < srf->w &&
					y < srf->h &&
					!((int)x == (int)xl)) {
					tosvg << x << "," << y << " ";
				}
				xl = x;
			}
			tosvg << "\" style=\"fill:#" + to_hexa(c.col) + ";fill-opacity:0.35\" />\n";
		}
	}
}

void Plotter::draw_main(const Curve & c) {

	double limright = srf->w - margin_right;
	double limdown = srf->h - margin_bottom;
	
	bool svg_safe = false;

	double xl = x_a * (c.x?c.x[0]:1) + x_b;
	double yl = y_a * c.y[0] + y_b;
	for(int i = 1; i < c.n; i++) {
		double x = x_a * (c.x?c.x[i]:(i+1)) + x_b;
		double y = y_a * c.y[i] + y_b;
		
		double xx1 = xl;
		double yy1 = yl;
		double xx2 = x;
		double yy2 = y;

		if(xx1 >= 0 &&
			yy1 >= 0 &&
			xx1 < srf->w &&
			yy1 < srf->h &&
			xx2 >= 0 &&
			yy2 >= 0 &&
			xx2 < srf->w &&
			yy2 < srf->h) {
			if(tosvg_on) {
				if(svg_safe) {
					if(!((int)xx1 == (int)xx2 /*&& (int)yy1 == (int)yy2*/))
						tosvg << " " << xx2 << "," << yy2;
				} else {
					tosvg << "<polyline points=\"" << xx1 << "," << yy1 << " " << xx2 << "," << yy2;
					svg_safe = true;
				}
			}
			draw_line(srf, (int)xx1, (int)yy1, (int)xx2, (int)yy2, c.col);
		}
		else {
			if(tosvg_on) {
				if(svg_safe) {
					tosvg << "\" style=\"fill:none;stroke:#"+to_hexa(c.col)+";stroke-width:1.5\" />\n";
					svg_safe = false;
				}
			}
		}

		xl = x;
		yl = y;
	}
	if(tosvg_on) {
		if(svg_safe) {
			tosvg << "\" style=\"fill:none;stroke:#"+to_hexa(c.col)+";stroke-width:1.5\" />\n";
			svg_safe = false;
		}
	}
		
}

void Plotter::drawVer(double x) {
	double xx = x_a * x + x_b;
	safe_draw_line(srf, xx, margin_top, xx, srf->h - margin_bottom, 0x808080);
}

void Plotter::drawHor(double y) {
	double yy = y_a * y + y_b;
	safe_draw_line(srf, margin_left, yy, srf->w - margin_right, yy, 0x808080);
}

//Replaces all occurrences of key in str by wrt
void str_replace(string & str, string key, string wrt) {
	while(true) {
		size_t found = str.find(key);
		if(found==string::npos)
			break;
		str.replace(str.find(key),key.length(),wrt);
	}
}

//Replaces occurrences of "_{###}" and "^{###}" (where "###" is some text) by the subscript and superscript markup codes of SVG
//Replaces "\\{" with "{" and "\\}" with "}"
string svg_process_text(string text) {
	str_replace(text, "_{", "<tspan style=\"font-size:65%;baseline-shift:sub\">");
	str_replace(text, "^{", "<tspan style=\"font-size:65%;baseline-shift:super\">");
	str_replace(text, "\\{", "{");
	str_replace(text, "}", "</tspan>");
	str_replace(text, "\\</tspan>", "}");
	return text;
}

void Plotter::drawText(int x, int y, const char * text, int mode, int fontsize) {
	SDL_Surface * surt = TTF_RenderText_Solid(font, text, (SDL_Color){0,0,0,0});
    SDL_Rect offset;
    offset.x = x;
    offset.y = y;
    bool rotate = false;
    switch(mode) {
    case 0: //origin at top-left
	    offset.x = x + 10;
	    offset.y = y + 2;
    	break;
    case 1: //origin at top-center
    	offset.x -= surt->w/2;
    	offset.y += 2;
    	break;
    case 2: //origin at center-right
    	offset.x -= surt->w + 2;
    	offset.y -= surt->h/2;
    	break;
    case 3: //origin at bottom-center, rotated
    	{
			SDL_Surface * new_surt = SDL_CreateRGBSurface(SDL_SWSURFACE, surt->h, surt->w, 32, 0, 0, 0, 0);
			for(int i = 0; i < surt->w; i++) {
				for(int j = 0; j < surt->h; j++) {
					((Uint32*)(new_surt->pixels))[j+(surt->w-1-i)*surt->h] = (1-((Uint8*)(surt->pixels))[i+j*(((surt->w+3)/4)*4)])*0xffffff;
				}
			}
			SDL_FreeSurface(surt);
			surt = new_surt;
    	}
    	offset.x -= surt->w + 4;
    	offset.y -= surt->h/2;
    	rotate = true;
    	break;
    case 4: //origin at center
    	offset.x -= surt->w/2;
    	offset.y -= surt->h/2;
    	break;
    }
	if(tosvg_on) {
		if(rotate) {
			std::string sx = to_stdstring(offset.x+surt->w);
			std::string sy = to_stdstring(offset.y+surt->h);
			tosvg << "  <text x=\"" + sx + "\" y=\"" + sy +"\" font-family=\"" + FONTNAME + "\" transform=\"rotate(-90 " + sx + "," + sy + ")\" font-size=\"" + to_stdstring(fontsize) + "\"> " + svg_process_text(text) + "</text>\n";
		}
		else
			tosvg << "  <text x=\""+to_stdstring(offset.x-1)+"\" y=\""+to_stdstring(offset.y+surt->h-3)+"\" font-family=\"" + FONTNAME + "\" font-size=\"" + to_stdstring(fontsize) + "\"> " + svg_process_text(text) + "</text>\n";
	}
    SDL_BlitSurface(surt, NULL, srf, &offset);
	SDL_FreeSurface(surt);
}

void Plotter::update_ab() {
	//COMPUTE MIN & MAX
	double xmin = xmin0;
	double ymin = ymin0;
	double xmax = xmax0;
	double ymax = ymax0;
	for(std::list<Curve>::iterator j = curves.begin(); j != curves.end(); j++) {
		Curve & c = *j;
		for(int i = 0; i < c.n; i++) {
			double x = c.x?c.x[i]:(i+1);
			double y = c.y[i];
			if(x-x == 0 && y-y==0) { //avoiding NaN
				if(x > xmax)
					xmax = x;
				if(x < xmin)
					xmin = x;
				if(y > ymax)
					ymax = y;
				if(y < ymin)
					ymin = y;
			}
		}
	}
	if(!(xmin - xmin == 0)) {
		x_a = y_a = x_b = y_b = 0;
		return;
	};
	
	//x_a * xmin + x_b = margin_left
	//x_a * xmax + x_b = srf->w - margin_right
	x_a = (srf->w - margin_right - margin_left) / (xmax-xmin); //positive
	y_a = (srf->h - margin_bottom - margin_top) / (ymin-ymax); //negative
	x_b = margin_left - x_a * xmin;
	y_b = margin_top - y_a * ymax;
	
	//TRUNCATE X
	double xsd_log = log(xmax - xmin) / log(10.);
	int xsd_log_floor = std::floor(xsd_log);
	double xsd_log_err = xsd_log - xsd_log_floor;
	int xsd_pref = 0;
	int xsd_pow = 0;
	if(xsd_log_err < log(2.)/log(10.)) {
		xsd_pow = xsd_log_floor - 1;
		xsd_pref = 2;
	}
	else if(xsd_log_err < log(5.)/log(10.)) {
		xsd_pow = xsd_log_floor - 1;
		xsd_pref = 5;
	}
	else {
		xsd_pow = xsd_log_floor;
		xsd_pref = 1;
	}
	double xsd = xsd_pref * pow(10, xsd_pow);
	//drawText(srf->w - margin_right, srf->h - margin_bottom, xaxisname.c_str(), 0, 21);
	drawText((margin_left + srf->w - margin_right)/2 + 30, srf->h - margin_bottom + 40, xaxisname.c_str(), 1, 25);

	//DRAW X LINES AND VALUES
	for(int xx = std::ceil(xmin / xsd); xx < xmax / xsd; xx++) {
		drawVer(xx * xsd);
		drawText(x_a * (xx*xsd) + x_b, srf->h - margin_bottom, to_stdstring(xsd * xx).c_str(), 1, 18);
	}
	
	//TRUNCATE Y
	double ysd_log = log(ymax - ymin) / log(10.);
	int ysd_log_floor = std::floor(ysd_log);
	double ysd_log_err = ysd_log - ysd_log_floor;
	int ysd_pref = 0;
	int ysd_pow = 0;
	if(ysd_log_err < log(2.)/log(10.)) {
		ysd_pow = ysd_log_floor - 1;
		ysd_pref = 2;
	}
	else if(ysd_log_err < log(5.)/log(10.)) {
		ysd_pow = ysd_log_floor - 1;
		ysd_pref = 5;
	}
	else {
		ysd_pow = ysd_log_floor;
		ysd_pref = 1;
	}
	double ysd = ysd_pref * pow(10, ysd_pow);
	//drawText(margin_left, margin_top, yaxisname.c_str(), 3, 21);
	drawText(margin_left - 60, (margin_top + srf->h - margin_bottom)/2 - 20, yaxisname.c_str(), 3, 25);

	//DRAW Y LINES AND VALUES
	for(int yy = std::ceil(ymin / ysd); yy < ymax / ysd; yy++) {
		drawHor(yy * ysd);
		drawText(margin_left, y_a * (yy*ysd) + y_b, to_stdstring(ysd * yy).c_str(), 2, 18);
	}
}

void Plotter::clearGraph() {
	for(int i = 0; i < srf->w*srf->h; i++)
		((Uint32*)(srf->pixels))[i] = 0xffffff;
	safe_draw_line(srf, margin_left, margin_top, margin_left, srf->h - margin_bottom, 0);
	safe_draw_line(srf, srf->w - margin_right, margin_top, srf->w - margin_right, srf->h - margin_bottom, 0);
	safe_draw_line(srf, margin_left, margin_top, srf->w - margin_right, margin_top, 0);
	safe_draw_line(srf, margin_left, srf->h - margin_bottom, srf->w - margin_right, srf->h - margin_bottom, 0);

	if(title != "")
		drawText((margin_left + srf->w - margin_right)/2, margin_top/2, title.c_str(), 4, 16);

}

void Plotter::drawEverything() {
	clearGraph();
	int counter = 0;
	for(std::list<Curve>::iterator j = curves.begin(); j != curves.end(); j++) {
		int hh = margin_top + 30 + (counter)*30;
		int hh2;
		bool tosvg_local = tosvg_on;
		tosvg_on = false;
		if(j->y_u && j->y_d) {
			hh2 = hh + 1;  safe_draw_line(srf, srf->w - margin_right + 20, hh2, srf->w - margin_right + 50, hh2, j->col_u);
			hh2 = hh + 11; safe_draw_line(srf, srf->w - margin_right + 20, hh2, srf->w - margin_right + 50, hh2, j->col_d);
			if(tosvg_local) {
				tosvg << "<polygon points=\"" << 
					srf->w - margin_right + 20 << "," << hh+1 << " " <<
					srf->w - margin_right + 50 << "," << hh+1 << " " <<
					srf->w - margin_right + 50 << "," << hh+11 << " " <<
					srf->w - margin_right + 20 << "," << hh+11 <<
				"\" style=\"fill:#" << to_hexa(j->col) << ";fill-opacity:0.35\" />\n";

			}
		}
		hh2 = hh + 6;  safe_draw_line(srf, srf->w - margin_right + 20, hh2, srf->w - margin_right + 50, hh2, j->col);
		if(tosvg_local) {
			tosvg << "<polyline points=\"" <<
				srf->w - margin_right + 20 << "," << hh+6 << " " <<
				srf->w - margin_right + 50 << "," << hh+6 <<
				"\" style=\"fill:none;stroke:#"+to_hexa(j->col)+";stroke-width:1.5\" />\n";
		}
		tosvg_on = tosvg_local;
		drawText(srf->w - margin_right + 51, hh-4, j->name.c_str(), 0, 18);
		counter++;
	}
	update_ab();
	for(std::list<Curve>::iterator j = curves.begin(); j != curves.end(); j++)
		draw_ud(*j);
	for(std::list<Curve>::iterator j = curves.begin(); j != curves.end(); j++)
		draw_main(*j);
}

void Plotter::refresh() {
	drawEverything();

	SDL_Event evt;
	while(SDL_PollEvent(&evt)) {
		switch(evt.type) {
		case SDL_QUIT:
			quit = true;
			break;
		case SDL_KEYDOWN:
			switch(evt.key.keysym.sym) {
			case SDLK_ESCAPE:
				quit = true;
				break;
			}
			break;
		}
	};
	
	if(quit)
		exportSVG();

	SDL_Flip(srf);
}

void Plotter::exportSVG() {
	cout << "Exporting " << svg_filename << "..." << endl;
	tosvg.str("");
	tosvg.clear();
	tosvg_on = true;
	tosvg << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
	tosvg << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
	drawEverything();
	tosvg << "</svg>\n";
	ofstream f(svg_filename.c_str());
	tosvg_on = false;
	f << tosvg.str();
	f.close();
}

Plotter::Curve * Plotter::addCurve() {
	curves.push_back(Curve());
	curves.back().name = "curve #" + to_stdstring(curves.size());
	return & curves.back();
};

Plotter::Plotter(int width, int height) {
	xmin0 = 1./0.;
	ymin0 = 1./0.;
	xmax0 = -1./0.;
	ymax0 = -1./0.;
	margin_left = 90;
	margin_right = 300;
	margin_top = 30;
	margin_bottom = 60;
	quit = false;
	svg_filename = "out.svg";
	xaxisname = "X";
	yaxisname = "Y";
	title = "";
	SDL_Init(SDL_INIT_EVERYTHING);
	srf = SDL_SetVideoMode(width, height, 32, SDL_SWSURFACE);
	TTF_Init();
	font = TTF_OpenFont(FONTPATH, 18);
	if(!font) {
		cout << "ERROR: COULD NOT LOAD FONT" << endl;
		exit(1);
	}
	tosvg_on = false;
}

Plotter::~Plotter() {
	TTF_CloseFont(font);
	TTF_Quit();
	SDL_FreeSurface(srf);
	SDL_Quit();
}

