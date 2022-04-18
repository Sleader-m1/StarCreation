#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QtCore>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    friend std::istream& operator>>(std::istream& stream, MainWindow& win);
    ~MainWindow();
    bool ellipse = false;
    const QString textB[2] = {"Start", "Stop"};
    QTimer *timer = new QTimer(this);
    QElapsedTimer Gtime;
private:
    Ui::MainWindow *ui;
protected:
    void paintEvent(QPaintEvent *event);

private slots:
    void on_bPulse_clicked();
    void on_MomentOfPulse_clicked();
    void on_bOutputStars_clicked();
    void on_bFile_clicked();
    void on_bEnergy_clicked();
};

#endif // MAINWINDOW_H
