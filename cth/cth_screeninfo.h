
#include <QApplication>
#include <QGuiApplication>
#include <QtWidgets/QSystemTrayIcon>
#include <QtWidgets/QDesktopWidget>
#include <QtWidgets/QMainWindow>

class cth_ScreenInfo: public QMainWindow
{
    Q_OBJECT
    
public:
    explicit cth_ScreenInfo(QWidget *parent = 0, QScreen* scr =0);
    virtual ~cth_ScreenInfo(); 

protected:
    bool event(QEvent *event) override;
    bool functionAfterShownCalled = false;

private:
    QScreen *curr_screen;

};
