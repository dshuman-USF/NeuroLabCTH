

#include "cthguidropedit.h"
#include <QFileInfo>

cthGuiDropEdit::cthGuiDropEdit (QWidget *parent) : QLineEdit(parent)
{
   setAcceptDrops(true);
}

cthGuiDropEdit::~cthGuiDropEdit()
{
}

// if here, we've decided to allow the drop
void cthGuiDropEdit::dropEvent(QDropEvent* event)
{
   QString text = event->mimeData()->text();
   text = text.remove("file://");
   QFileInfo info(text);
   QString name = info.fileName();
   setText(name);
   event->acceptProposedAction();
}

// user proposed droping something in our lineedit. We we allow?
void cthGuiDropEdit::dragEnterEvent(QDragEnterEvent *event)
{
    if (event->mimeData()->hasText() ||
        event->mimeData()->hasUrls())
    {
      QString text = event->mimeData()->text();
      text = text.remove("file://"); // turn url into path/fname
      QFileInfo info(text);
      if (info.isFile() && info.completeSuffix() == "type")
         event->acceptProposedAction();
    }
}

