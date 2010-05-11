/****************************************************************************
**
** Copyright (C) 2005-2005 Trolltech AS. All rights reserved.
**
** This file is part of the example classes of the Qt Toolkit.
**
** This file may be used under the terms of the GNU General Public
** License version 2.0 as published by the Free Software Foundation
** and appearing in the file LICENSE.GPL included in the packaging of
** this file.  Please review the following information to ensure GNU
** General Public Licensing requirements will be met:
** http://www.trolltech.com/products/qt/opensource.html
**
** If you are unsure which license is appropriate for your use, please
** review the following information:
** http://www.trolltech.com/products/qt/licensing.html or contact the
** sales department at sales@trolltech.com.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
****************************************************************************/

/*
    delegate.cpp

    A delegate that allows the user to change integer values from the model
    using a spin box widget.
*/

#include <QtGui>

#include "delegate.h"


SpinBoxDelegate::SpinBoxDelegate(QObject *parent)
    : QItemDelegate(parent)
{
	currentEditor = NULL;
	bIsEditorCreated = false;
}

QWidget *SpinBoxDelegate::createEditor(QWidget *parent,
    const QStyleOptionViewItem &/* option */,
    const QModelIndex & index )  const 
{
    QSpinBox *editor = new QSpinBox(parent);
    editor->setMinimum(0);
    editor->setMaximum(100);
    editor->installEventFilter(const_cast<SpinBoxDelegate*>(this));

    
    const_cast<SpinBoxDelegate*>(this)->setCurrentEditor( editor, index );

    return editor;
}

void SpinBoxDelegate::setCurrentEditor(  const QWidget * editor, const QModelIndex index )
{
	bIsEditorCreated = true;
	currentEditor = const_cast< QWidget * >(editor);
//	currentModelIndex = const_cast<QModelIndex>(index);
	currentModelIndex = index;

}

void SpinBoxDelegate::setEditorData(QWidget *editor,
                                    const QModelIndex &index) const
{
    int value = index.model()->data(index, Qt::DisplayRole).toInt();

    QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
    spinBox->setValue(value);
}

void SpinBoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
    QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
    spinBox->interpretText();
    int value = spinBox->value();

    model->setData(index, value);
    const_cast<SpinBoxDelegate*>(this)->setIsEditorCreated( false );
}

void SpinBoxDelegate::updateEditorGeometry(QWidget *editor,
    const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
    editor->setGeometry(option.rect);
}

void SpinBoxDelegate::commitCurrentEditor( QAbstractItemModel * model)
{
	if( currentEditor == NULL || !bIsEditorCreated )
		return;

	setModelData( currentEditor, model, currentModelIndex );
	
}

void SpinBoxDelegate::setIsEditorCreated( bool bEditor )
{
	bIsEditorCreated = bEditor;

}
