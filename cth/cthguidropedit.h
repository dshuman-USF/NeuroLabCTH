#ifndef CTHGUIDORPEDIT_H
#define CTHGUIDORPEDIT_H

#include <QObject>
#include <QLineEdit>
#include <QDropEvent>
#include <QDragEnterEvent>
#include <QMimeData>
#include <QValidator>

class cthGuiDropEdit : public QLineEdit
{
    Q_OBJECT
    
   public:
    explicit cthGuiDropEdit(QWidget *parent = 0);
    virtual ~cthGuiDropEdit();

   protected:
    void dropEvent(QDropEvent*) override;
    void dragEnterEvent(QDragEnterEvent *event) override;
};

#endif
