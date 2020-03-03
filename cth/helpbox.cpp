
#include <QDesktopServices>
#include "helpbox.h"
#include "ui_helpbox.h"
#include <iostream>

helpbox::helpbox(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::helpbox)
{
    ui->setupUi(this);
}

helpbox::~helpbox()
{
    delete ui;
}


void helpbox::on_helpClose_accepted()
{
   close();
}

void helpbox::on_helpClose_rejected()
{
   close();
}

void helpbox::on_helpText_anchorClicked(const QUrl &arg1)
{
   QDesktopServices::openUrl(QUrl(arg1));
}
