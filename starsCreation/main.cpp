#include "mainwindow.h"
#include <QApplication>
#include <sstream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    //std::string str = "1000 1000 1 -1.28581e+11 -1.28581e+11 10 10 5.9742e24";//0 100 100 100 100 10000000 23 24 100 100 100000 60 70 100 100 100000000 10 23 100 100 120000 14 45 100 100 320000 34 53 100 100 10000 34 53 100 100 400000 75 35 100 100 10000 34 54 100 100 10000 13 53 100 100 200000";
    //std::stringstream ss (str);
    //ss >> w;
    w.show();
    return a.exec();
}
