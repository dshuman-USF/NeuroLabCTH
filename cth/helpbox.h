#ifndef HELPBOX_H
#define HELPBOX_H

#include <QDialog>

namespace Ui {
class helpbox;
}

class helpbox : public QDialog
{
    Q_OBJECT

public:
    explicit helpbox(QWidget *parent = 0);
    virtual ~helpbox();

private slots:
    void on_helpClose_accepted();
    void on_helpClose_rejected();
    void on_helpText_anchorClicked(const QUrl &arg1);

private:
    Ui::helpbox *ui;
};

#endif // HELPBOX_H
