#include "cmath"
#include <QPainter>
#include <QFile>
#include <iomanip>
#include <QTextStream>
#include <sstream>

#ifndef STAR_H
#define STAR_H
const int dim = 2;
int numStars = 100;
double systemRadius = 1e12, distConnect = 1e11;
const int borderMassC = 10;
const double G = 6.67408e-11, //расстояние для слияния
dt = 10000;
const double massSun   = 1.98892e30,
             massJup   = 1898.6e24,
             massUran  = 86.832e24,
             massEarth = 5.9742e24,
             massVenus = 4.867e24;
const double borderMass[] = {borderMassC*massEarth, borderMassC*massUran, borderMassC*massJup, borderMassC*massSun};
const QColor colStar[] = {Qt::cyan, Qt::darkGreen, Qt::magenta, Qt::yellow, Qt::white};
const int nColor = sizeof(colStar) / sizeof(colStar[0]);

class star{
public:
    static int starCounter;
    double x[dim];
    double v[dim];
    double m;
    double f[dim];
    QColor col;
    star(double *coord, double *speed, double mass);
    friend std::ostream& operator<<(std::ostream& stream, star& s);
    friend star operator+(star str1, star str2);
    ~star(){starCounter--;}
};
int star::starCounter = 0;

std::ostream& operator<<(std::ostream& stream, star& s){
    stream << "(";
    for(int i = 0; i < dim; ++i){
        stream << std::setw(15) << std::setprecision(8) <<std::setfill(' ') << std::right  << s.x[i];
        if(i != dim -1) stream << "; ";
    }
    stream << ") ";
    stream << "vx = " << std::setw(8)<< std::setprecision(6) << std::right <<  s.v[0] << ", ";
    stream << "vy = "<< std::setw(8)<< std::setprecision(6) << std::right << s.v[1];
    if(dim == 3){
        stream << ", vz = " << s.v[2];
    }
    stream << ", ";
    stream << "mass = " << s.m;
    return stream;
}

std::istream& operator>>(std::istream& stream, star& s){
    stream >> s.x[0];
    stream >> s.x[1];
    stream >> s.v[0];
    stream >> s.v[1];
    stream >> s.m;
    return stream;
}

star operator+(star str1, star str2){
    for(int i = 0; i < dim; ++i){
        str1.v[i] = (str1.v[i]*str1.m + str2.v[i]*str2.m)/(str1.m+str2.m);
        str1.x[i] = (str1.x[i]*str1.m + str2.x[i]*str2.m)/(str1.m+str2.m);
    }
    str1.m += str2.m;
    return str1;
}

star::star(double *coord, double *speed, double mass){
    for(int k = 0; k < dim; ++k){
        x[k] = coord[k];
        v[k] = speed[k];
    }
    m = mass;
    col = colStar[nColor-1];  // не годится, если будут объединения объектов, функция намного лучше
    for(int i = 0; i < nColor-1; ++i){
        if(m <= borderMass[i]){
            col = colStar[i];
            break;
        }
    }
    starCounter++;
}
class galaxy{
    long double massOfAllSystem;
public:
    int num;
    star **stars;

    double getDist(double* x1, double* x2){
        return sqrt((x1[0]+x2[0])*(x1[0]+x2[0])+(x1[1]+x2[1])*(x1[1]+x2[1]));
    }

    QString GetStarsCount(){
        int counter = 0;
        for(int i = 0; i < num; ++i){
            if (stars[i] != nullptr)
                ++counter;
        }
        std::stringstream ss;
        ss << counter;
        return QString::fromStdString(ss.str());
    }

    void setGalaxy(int sNum, double sRad, double dConnect){
        numStars = sNum;
        systemRadius = sRad;
        distConnect = dConnect;
    }

    void SumStars(){
        for(int i = 1; i < num-1; ++i){
            if(stars[i] != nullptr)
                for(int j = i+1; j < num; ++j){
                    if(stars[j] != nullptr){
                        if(abs(stars[i]->x[0] - stars[j]->x[0]) < distConnect){
                            if(abs(stars[i]->x[1] - stars[j]->x[1]) < distConnect){
                                *stars[i] = *stars[i] + *stars[j];
                                stars[j] = nullptr;
                            }
                        }
                    }
                }
        }
    }

    QString GetMassOfAllGalaxy(){
        std::stringstream ss;
         ss << massOfAllSystem;
        return QString::fromStdString(ss.str());
    }

    QString GetMomentPulseX(){
        long double MpulseX = 0;
        for(int i = 1; i < num; ++i){
                if(stars[i] != nullptr){
                    double dist = std::sqrt((stars[0]->x[0] - stars[i]->x[0])*(stars[0]->x[0] - stars[i]->x[0])+(stars[0]->x[1] - stars[i]->x[1])*(stars[0]->x[1] - stars[i]->x[1]));
                MpulseX += stars[i]->v[0]*(stars[i]->m/massSun)*dist;
                }
           }
        std::stringstream ss;
         ss << MpulseX;
        return QString::fromStdString(ss.str());
    }

    QString GetMomentPulseY(){
        long double MpulseY = 0;
        for(int i = 1; i < num; ++i){
            if(stars[i] != nullptr){
                double dist = std::sqrt((stars[0]->x[0] - stars[i]->x[0])*(stars[0]->x[0] - stars[i]->x[0])+(stars[0]->x[1] - stars[i]->x[1])*(stars[0]->x[1] - stars[i]->x[1]));
                MpulseY += stars[i]->v[1]*(stars[i]->m/massSun)*dist;
            }
        }
        std::stringstream ss;
         ss << MpulseY;
        return QString::fromStdString(ss.str());
    }

    QString GetPulseX(){
        long double pulseX = 0;
        for(int i = 1; i < num; ++i){
                if(stars[i] != nullptr)
                pulseX += stars[i]->v[0]*(stars[i]->m/massSun);
           }
        std::stringstream ss;
         ss << pulseX;
        return QString::fromStdString(ss.str());
    }

    QString GetPulseY(){
        long double pulseY = 0;
        for(int i = 1; i < num; ++i){
            if(stars[i] != nullptr)
                pulseY += stars[i]->v[1]*(stars[i]->m/massSun);
        }
        std::stringstream ss;
         ss << pulseY;
        return QString::fromStdString(ss.str());
    }

    double GetEnergyX(){
        long double energyX = 0;
        for(int i = 1; i < num; ++i){
                if(stars[i] != nullptr)
                energyX += stars[i]->v[0]*stars[i]->v[0]*stars[i]->m/2;
           }
        return energyX;
    }

    double GetEnergyY(){
        long double energyY = 0;
        for(int i = 1; i < num; ++i){
            if(stars[i] != nullptr)
                energyY += stars[i]->v[1]*stars[i]->v[1]*stars[i]->m/2;
        }
        return energyY;
    }

    galaxy(int n = numStars):num(n){
        stars = new star*[num];
        double x1[dim] = {0}, v1[dim] = {0};
        stars[0] = new star(x1, v1, massSun); // самый массивный объект в начале координат
        massOfAllSystem = massSun;
        double rad;
        for(int i = 1; i < num; ++i){
            rad = 0;
            double R = rand() * systemRadius / RAND_MAX,
            fi  = (2 * M_PI * rand()) / RAND_MAX,
            theta = (M_PI * rand()) / RAND_MAX;
            x1[0] = R  * cos(fi);
            x1[1] = R  * sin(fi);
            if(dim == 3){
                x1[0] *= sin(theta);
                x1[1] *= sin(theta);
                x1[3] = R * cos(theta);
            }
            for(int k = 0; k < dim; ++k){
                rad += x1[k] * x1[k];
            }
            // вторая космическая скорость
            double absV = sqrt(G * stars[0]->m / sqrt(rad)), alpha = (2 * M_PI * rand()) / RAND_MAX;

            v1[0] =  absV * sin(fi);
            v1[1] = -absV * cos(fi); // скорость направлена вдоль окружности с центром в начале координат
            stars[i] = new star(x1, v1, massEarth / num * (6 * i));
            massOfAllSystem += stars[i]->m;
        }
    };

    friend std::ostream& operator<<(std::ostream& stream, galaxy& galac){
        for(int i = 0; i < galac.num; ++i){
            if (galac.stars[i] != nullptr)
            stream << "star["<<std::right<< std::setw(3) << std::setfill('0')<< i << "]" << *galac.stars[i] << std::endl;
        }
        return stream;
    }

    friend std::istream& operator>>(std::istream& stream, galaxy& galac){
        stream >> galac.num;
        galac.num++;
        double x1[dim] = {0};
        double v1[dim] = {0};
        galac.stars[0] = new star(x1, v1, massSun);
        for(int i = 1; i < galac.num; ++i){
            stream >> *galac.stars[i];
        }
        return stream;
    }

    ~galaxy(){delete[] stars;};
    double border[dim];

    void move(){
        double dist;
        double dCoord[dim];
        for(int i = 0; i < num; ++i){ // force nullification
            for(int k = 0; k < dim; ++k){
                if(stars[i]){
                    stars[i]->f[k] = 0;
                }
            }
        }
        for(int i = 0; i < num; i++){
            if(stars[i]){
                for(int j = i + 1; j < num; j++){
                    if(i != j && stars[j]){
                        dist = 0;
                        for(int k = 0; k < dim; ++k){
                            dCoord[k] = stars[i]->x[k] - stars[j]->x[k];
                            dist += dCoord[k] * dCoord[k];
                        }
                        if(sqrt(dist) < distConnect){
                            double tmpM = stars[i]->m + stars[j]->m, tmpX[dim], tmpV[dim];
                            for(int k = 0; k < dim; ++k){
                                tmpX[k] = (stars[i]->x[k] * stars[i]->m + stars[j]->x[k] * stars[j]->m)/tmpM;
                                tmpV[k] = (stars[i]->v[k] * stars[i]->m + stars[j]->v[k] * stars[j]->m)/tmpM;
                            }
                            delete stars[j];
                            stars[j] = nullptr;
                            stars[i]->m = tmpM;
                            for(int k = 0; k < dim; ++k){
                                stars[i]->x[k] = tmpX[k];
                                stars[i]->v[k] = tmpV[k];
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num; i++){
            if(stars[i]){
                for(int j = i + 1; j < num; j++){
                    if(i != j && stars[j]){
                        dist = 0;
                        for(int k = 0; k < dim; ++k){
                            dCoord[k] = stars[i]->x[k] - stars[j]->x[k];
                            dist += dCoord[k] * dCoord[k];
                        }
                        // dist = sqrt(dist); // для знаменателя пока квадрат, потом возьмем корень
                        for(int k = 0; k < dim; ++k){
                            double tmp = G * stars[i]->m * stars[j]->m / dist;
                            stars[i]->f[k] -= tmp * dCoord[k] / sqrt(dist);
                            stars[j]->f[k] += tmp * dCoord[k] / sqrt(dist);
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num; ++i){
            if(stars[i]){
                for(int k = 0; k < dim; ++k){
                    stars[i]->v[k] += dt * stars[i]->f[k] / stars[i]->m; //можно не делить на массу, а выше суммировать ускорение
                }
                for(int k = 0; k < dim; ++k){
                    stars[i]->x[k] += dt * stars[i]->v[k];
                }
            }
        }
    }
};

#endif // STAR_H
