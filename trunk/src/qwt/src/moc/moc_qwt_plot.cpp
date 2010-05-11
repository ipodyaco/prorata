/****************************************************************************
** Meta object code from reading C++ file 'qwt_plot.h'
**
** Created: ??? ?? 25 15:38:30 2006
**      by: The Qt Meta Object Compiler version 59 (Qt 4.1.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../include/qwt_plot.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qwt_plot.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.1.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_QwtPlot[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   10, // methods
       7,   45, // properties
       1,   66, // enums/sets

 // signals: signature, parameters, type, tag, flags
      18,    9,    8,    8, 0x05,
      58,   46,    8,    8, 0x05,

 // slots: signature, parameters, type, tag, flags
      91,    8,    8,    8, 0x0a,
      99,    8,    8,    8, 0x0a,
     108,    8,    8,    8, 0x0a,
     122,    8,    8,    8, 0x09,
     142,    8,    8,    8, 0x09,

 // properties: name, type, flags
     171,  166, 0x01095103,
     189,  182, 0x43095103,
     210,  206, 0x02095103,
     226,  166, 0x01095003,
     238,  166, 0x01095003,
     247,  166, 0x01095003,
     257,  166, 0x01095003,

 // enums: name, flags, count, data
     268, 0x0,    5,   70,

 // enum data: key, value
     273, uint(QwtPlot::yLeft),
     279, uint(QwtPlot::yRight),
     286, uint(QwtPlot::xBottom),
     294, uint(QwtPlot::xTop),
     299, uint(QwtPlot::axisCnt),

       0        // eod
};

static const char qt_meta_stringdata_QwtPlot[] = {
    "QwtPlot\0\0plotItem\0legendClicked(QwtPlotItem*)\0plotItem,on\0"
    "legendChecked(QwtPlotItem*,bool)\0clear()\0replot()\0autoRefresh()\0"
    "legendItemClicked()\0legendItemChecked(bool)\0bool\0autoReplot\0QColor\0"
    "canvasBackground\0int\0canvasLineWidth\0xBottomAxis\0xTopAxis\0"
    "yLeftAxis\0yRightAxis\0Axis\0yLeft\0yRight\0xBottom\0xTop\0axisCnt\0"
};

const QMetaObject QwtPlot::staticMetaObject = {
    { &QFrame::staticMetaObject, qt_meta_stringdata_QwtPlot,
      qt_meta_data_QwtPlot, 0 }
};

const QMetaObject *QwtPlot::metaObject() const
{
    return &staticMetaObject;
}

void *QwtPlot::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_QwtPlot))
	return static_cast<void*>(const_cast<QwtPlot*>(this));
    if (!strcmp(_clname, "QwtPlotDict"))
	return static_cast<QwtPlotDict*>(const_cast<QwtPlot*>(this));
    return QFrame::qt_metacast(_clname);
}

int QwtPlot::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QFrame::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: legendClicked((*reinterpret_cast< QwtPlotItem*(*)>(_a[1]))); break;
        case 1: legendChecked((*reinterpret_cast< QwtPlotItem*(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 2: clear(); break;
        case 3: replot(); break;
        case 4: autoRefresh(); break;
        case 5: legendItemClicked(); break;
        case 6: legendItemChecked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        }
        _id -= 7;
    }
#ifndef QT_NO_PROPERTIES
      else if (_c == QMetaObject::ReadProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< bool*>(_v) = autoReplot(); break;
        case 1: *reinterpret_cast< QColor*>(_v) = canvasBackground(); break;
        case 2: *reinterpret_cast< int*>(_v) = canvasLineWidth(); break;
        case 3: *reinterpret_cast< bool*>(_v) = xBottomAxisEnabled(); break;
        case 4: *reinterpret_cast< bool*>(_v) = xTopAxisEnabled(); break;
        case 5: *reinterpret_cast< bool*>(_v) = yLeftAxisEnabled(); break;
        case 6: *reinterpret_cast< bool*>(_v) = yRightAxisEnabled(); break;
        }
        _id -= 7;
    } else if (_c == QMetaObject::WriteProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: setAutoReplot(*reinterpret_cast< bool*>(_v)); break;
        case 1: setCanvasBackground(*reinterpret_cast< QColor*>(_v)); break;
        case 2: setCanvasLineWidth(*reinterpret_cast< int*>(_v)); break;
        case 3: enableXBottomAxis(*reinterpret_cast< bool*>(_v)); break;
        case 4: enableXTopAxis(*reinterpret_cast< bool*>(_v)); break;
        case 5: enableYLeftAxis(*reinterpret_cast< bool*>(_v)); break;
        case 6: enableYRightAxis(*reinterpret_cast< bool*>(_v)); break;
        }
        _id -= 7;
    } else if (_c == QMetaObject::ResetProperty) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 7;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void QwtPlot::legendClicked(QwtPlotItem * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void QwtPlot::legendChecked(QwtPlotItem * _t1, bool _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
