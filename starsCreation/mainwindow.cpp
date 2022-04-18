#include "mainwindow.h"
#include "math.h"

#include "ui_mainwindow.h"
#include <QPainter>
#include <QTime>
#include <star.h>
#include <fstream>

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow){
    ui->setupUi(this);
    ui->bPulse->click();
    Gtime.start();
    ui->MomentOfPulse->click();
    ui->bEnergy->click();
//    connect(ui->pushButtonStart,  SIGNAL(clicked()), this, SLOT(buttonText()));
    connect(timer, &QTimer::timeout, this, QOverload<>::of(&MainWindow::update));
//    connect(timer, &QTimer::timeout, this, SLOT(myTimer()));
    timer->start(1);
}
const int topX0 = 100, topY0 = 100;
int h = 800, length = 900;
galaxy *galactika = new galaxy;

std::istream& operator>>(std::istream& stream, MainWindow& win){
    int starsNum;
    double sysRad, distConn;
    stream >> starsNum >> sysRad>>distConn ;
    galactika->setGalaxy(starsNum, sysRad, distConn);
    galactika = new galaxy;
    return stream;
}

MainWindow::~MainWindow(){
    delete ui;
}
void MainWindow::paintEvent(QPaintEvent *e) {
  Q_UNUSED(e);
  QPainter painter(this);
  QPen pen(Qt::black, 1, Qt::SolidLine);
  painter.setPen(pen);
  QTime time = QTime::currentTime();
  QBrush brush;//(Qt::yellow);
  brush.setStyle(Qt::SolidPattern);
  double coefX = length / 2 / 1e12;
  int centerX = length / 2;
  for(int i = 0; i < galactika->num; ++i){
      if(galactika->stars[i]){
          brush.setColor(galactika->stars[i]->col);
          if(!i) brush.setColor(Qt::yellow);
          painter.setBrush(brush);
          for(int k = 0; k < 2; ++k){
                    painter.drawEllipse(galactika->stars[i]->x[0] * coefX + centerX + topX0,
                                        galactika->stars[i]->x[1] * coefX + centerX + topY0,
                                        6 + 4 * !i, 6 + 4 * !i);
          }
      }
  }
  galactika->move();
  if(time.second() % 5 == 0)
  galactika->SumStars();
  ui->leTimeField->setText("Current time: "+time.toString("hh:mm:ss:zzz"));
  QTime Htime(0, 0, 0);
  ui->bGalaxyTime->setText("Galaxy time: "+Htime.addMSecs(Gtime.elapsed()).toString("hh:mm:ss:zzz"));
  ui->leMassOfSystem->setText("Mass of galaxy = "+galactika->GetMassOfAllGalaxy());
  ui->bCounter->setText(galactika->GetStarsCount());
}

void MainWindow::on_bPulse_clicked()
{
    ui->lePulseX->setText("Px / massSun = "+galactika->GetPulseX());
    ui->lePulseY->setText("Py / massSun = "+galactika->GetPulseY());
}


void MainWindow::on_MomentOfPulse_clicked()
{
    ui->leMomentOfPulseX->setText("Lx = "+galactika->GetMomentPulseX());
    ui->leMomentOfPulseY->setText("Ly = "+galactika->GetMomentPulseY());
}

void MainWindow::on_bOutputStars_clicked()
{
    std::ofstream file("C:\\galaxy\\some.txt"); //Введите название файла для вывода
    file << *galactika << std::endl;
}


void MainWindow::on_bFile_clicked()
{
    std::ifstream file("C:\\galaxy\\some.txt");
    file >> *this;
    file.close();
}



void MainWindow::on_bEnergy_clicked()
{
    ui->leEnergyX->setText("Ex = " + QString::number(galactika->GetEnergyX()));
    ui->leEnergyY->setText("Ey = " + QString::number(galactika->GetEnergyY()));
}

