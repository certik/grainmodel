#ifndef SERIES_H
#define SERIES_H

class Updater
{
    public:
        virtual ~Updater() {}
        virtual void init(int maxval)=0;
        virtual void update(int i)=0;
};

void mesh(const std::string& fmesh, const std::string& fmatrices,
    const std::string& fboundaries, Updater *up);

void grad(const std::string& meshfile, double* x, int xsize, 
        double* g, int gsize, Updater *up);

double integ(const std::string& meshfile, const std::string& fboundaries, 
        double* x, int xsize, int b, Updater *up);

class loadmatrices
//currently Ax=F
//supports loading A and F from a file
{
    public:
        loadmatrices(const std::string& fname);
        ~loadmatrices();
        double readfloat();
        unsigned int readint();
    private:
        std::ifstream *f;
};


#endif
