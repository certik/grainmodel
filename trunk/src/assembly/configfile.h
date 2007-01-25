/* Quantum */
/* Copyright (C) 2006  Ondrej Certik */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include "tinyxml.h"

enum pottype {POT_WELL, POT_LHO, POT_HYDROGEN, POT_DFT};

class Config
{
public:
    Config(const char *fname)
    {
        nodesmap=true;
        eigs=3;
        solve=true;
        assemble=false;
        potential=POT_WELL;
        hbar=1.0; //atomic units
        m=1.0; //atomic units
        lhofrequency=1.0;
        dim=2;
        printlog=false;
        Eshift=0.0;
        scale=1.0;

        doc=new TiXmlDocument(fname);
    }
    ~Config()
    {
        delete doc;
    }
    int load()
    {
        doc->LoadFile();
	    TiXmlElement* rootElement=doc->FirstChildElement("root");
        if (!rootElement) 
        {
            std::cout << "Configuration file " << doc->Value() << " not found."
                << std::endl;
            error();
        }
        parse(rootElement);
        return 0;
    }
    void parseconfig(TiXmlElement *el)
    {
        TiXmlNode* c;
        c=el->FirstChild();
        while (c!=NULL)
        {
            const char* s=c->FirstChild()->Value();
            if (!strcmp(s,"print nodes map"))
                nodesmap=true;
            else if (!strcmp(s,"don't solve"))
                solve=false;
            else if (!strcmp(s,"assemble"))
                assemble=true;
            else if (!strcmp(s,"print log"))
                printlog=true;
            else
                std::cout << "warning - unknown option: " << s << std::endl;
            c=el->IterateChildren(c);
        }
    }
    void parse(TiXmlElement *el)
    {
        TiXmlNode* c;
        c=el->FirstChild();
        while (c!=NULL)
        {
            const char* s=c->Value();
            if (!strcmp(s,"config"))
                parseconfig(c->ToElement());
            else if (!strcmp(s,"eigenvalues"))
                eigs=atoi(c->FirstChild()->Value());
            else if (!strcmp(s,"dim"))
                dim=atoi(c->FirstChild()->Value());
            else if (!strcmp(s,"hbar"))
                hbar=atof(c->FirstChild()->Value());
            else if (!strcmp(s,"m"))
                m=atof(c->FirstChild()->Value());
            else if (!strcmp(s,"eshift"))
                Eshift=atof(c->FirstChild()->Value());
            else if (!strcmp(s,"scale"))
                scale=atof(c->FirstChild()->Value());
            else if (!strcmp(s,"lhofrequency"))
                lhofrequency=atof(c->FirstChild()->Value());
            else if (!strcmp(s,"potential"))
            {
                const char* s2=c->FirstChild()->Value();
                if (!strcmp(s2,"well"))
                    potential=POT_WELL;
                else if (!strcmp(s2,"lho"))
                    potential=POT_LHO;
                else if (!strcmp(s2,"hydrogen"))
                    potential=POT_HYDROGEN;
                else if (!strcmp(s2,"dft"))
                    potential=POT_DFT;
                else
                    std::cout << "warning - unknown option: " << 
                        s << ": " << s2 << std::endl;
            }
            else
                std::cout << "warning - unknown option: " << s << std::endl;
            c=el->IterateChildren(c);
        }
    }
    void print()
    {
        std::cout << "dimension: " << dim << std::endl;
        std::cout << "eigs: " << eigs << std::endl;
        std::cout << "nodesmap: " << nodesmap << std::endl;
        std::cout << "solve: " << solve << std::endl;
        std::cout << "assemble: " << assemble << std::endl;
        std::cout << "printlog: " << printlog << std::endl;
        std::cout << "potential: "; 
        switch (potential) {
            case POT_WELL: std::cout << "well" << std::endl; break;
            case POT_LHO: std::cout << "lho" << std::endl; break;
            case POT_HYDROGEN: std::cout << "hydrogen" << std::endl; break;
            default: std::cout << potential << std::endl; break;
        }
        if (potential==POT_LHO)
            std::cout << "lhofrequency: " << lhofrequency << std::endl;
        std::cout << "hbar: " << hbar << std::endl;
        std::cout << "m: " << m << std::endl;
        std::cout << std::endl;
    }

    int eigs;
    int dim;
    bool nodesmap;
    bool solve;
    bool assemble;
    bool printlog;
    pottype potential;
    double hbar,m;
    double Eshift,scale;
    double lhofrequency;
private:
    TiXmlDocument* doc;
};
