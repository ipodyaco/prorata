/****************************************************************************
** Meta object code from reading C++ file 'qwt_scale_widget.h'
**
** Created: ??? ?? 25 15:38:41 2006
**      by: The Qt Meta Object Compiler version 59 (Qt 4.1.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../include/qwt_scale_widget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qwt_scale_widget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.1.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_QwtScaleWidget[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_QwtScaleWidget[] = {
    "QwtScaleWidget\0"
};

const QMetaObject QwtScaleWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_QwtScaleWidget,
      qt_meta_data_QwtScaleWidget, 0 }
};

const QMetaObject *QwtScaleWidget::metaObject() const
{
    return &staticMetaObject;
}

void *QwtScaleWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_QwtScaleWidget))
	return static_cast<void*>(const_cast<QwtScaleWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int QwtScaleWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
