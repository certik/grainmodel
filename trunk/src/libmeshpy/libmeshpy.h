#ifndef SERIES_H
#define SERIES_H

double doubleSum( double* series, int size);
void doubleOnes(  double* array, int size);

void grad(const std::string& meshfile, double* x, int xsize, 
        double* g, int gsize);

void mesh(const std::string& meshfile);

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
