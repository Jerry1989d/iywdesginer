# editorpro

## 项目安装
```
npm install
```

### 开启本地服务
```
npm run serve
```

### 项目打包命令
```
npm run build
```

## Vue环境安装

vue 环境安装只需要引入本Demo直接输如命令安装即可
```jvascript
 npm i
```

#### 依赖包：
如需裁剪功能需要引入 [cropperjs](https://github.com/fengyuanchen/cropperjs)

3D功能需要引入 [Three.js](https://github.com/mrdoob/three.js/)

设计器插件单独安装

````javascript

npm i iywdesigner

````

引入用于合成的插件 
```html
<script src="https://cdn.bootcdn.net/ajax/libs/dom-to-image/2.6.0/dom-to-image.min.js"></script>
```


引入设计器插件

```javascript
import IYWDesigner from "./IYWDesigner";
```

### 设计器的调用
1.在需要生成设计器的组件里首先要给一个父级容器，并且设置一个id
2.设计器的默认定位是绝对定位，所以建议在editor-box外面添加一层相对定位的父级容器
```html
<div class="editor-box" id="editor"></div>
```
2.如需要预览效果则需要同上一样添加一个预览区的容器
```html
  <div class="preview-box" id="previewBox" draggable="false">
    <span class="preview-tag">预览</span>
</div>
```
在dom完成渲染后生成设计器
```javascript
    this.editorObj = new IYWDesigner(
        {
          id:'editor',// 父级容器的ID
          previewBox:'previewBox', // 预览区域的ID
          selectedInfo:_this.selectedInfo, // 如果需要得到当前操作对象的信息则传
          exportSize: 800, // 导出图的尺寸
        },
    );
    this.editorObj.init(); // 初始化设计器
```

## Options
初始化插件可以允许你传入一些选项值来自定义设计器的功能。

### id
  * Type: ``String`` * 必传项

### previewBox
* Type: ``String`` * 必传项

### selectedInfo
* Type: ``Object``  该对象下必须包含一个``info`` 对象用来接受设计器内当前选中对象的属性
* Return: 返回当前选中对象的以下属性 ``{width, height, zIndex, type, left, right}`` 


### maxItemNum
* Type: ``Number``' 编辑器内包含图片和色块的最大数量，默认为6 

### colorMax
* Type: ``Number``' 色块的最大数量限制，默认为5

## Methods


### init()
在页面生成设计器


```javascript
const editor = new IYWDesigner(options);
editor.init()
```
### addImage(imageInfo)
添加图案到设计器
```javascript
const imageInfo = { name: '示例图', type: 'pattern', id: 1, src:'1.png' };
editor.addImage(imageInfo)
```
* name: ``String`` 图片的名称 
* type: ``String`` 图片类型 1: 授权. 2: 素材. 3: IP. 4: 本地 (如果此项没有传默认识别为本地图案)
* id: ``Number`` 图案的id
* src: ``String`` * **图片的url** （必传项）
### delElement()
删除当前选中的对象
```javascript
editor.delElement()
```

### setTemp(sampleInfo)
设置样版
```javascript
const sampleInfo = {
    name: '示例样板',
    id: 1, 
    src: "../assets/static/images/sample-1.png",
    bladeInfo:[
      {
        points: [],  
        name:"blade1",
        bladeImg: "https://chdesign.oss-cn-shanghai.aliyuncs.com/test/CUMS/OA/20210225/161425967385244942.png",
      }
    ] 
};
editor.blad.set(0); // 注意这里一定要设置
editor.setTemp(sampleInfo);
```
* name: ``String`` 样版的名称
* id: ``Number`` 样版的id
* src: ``String`` * **样版的url** （必传项）
* bladeInfo: ``Array`` 
  * name : ``String`` 刀片的名称
  * bladeImg: ``String`` 刀片图地址
  * points: ``Array`` 该刀片图贴图区域的坐标点
  

### toolKit(key)
设计器操作工具集合函数，接受一个参数
```javascript
editor.toolKit('layerUp');
```

#### key值对照表
  * layerUp: 上一层
  * layerDown: 下一层
  * centerW: 水平居中
  * centerH: 垂直居中
  * mirrorX: 水平镜像
  * mirrorY: 垂直镜像
  * maxWidth: 宽度撑满
  * maxHeight: 高度撑满
  * fullCover: 铺满画布
  * clip: 对图像进行裁剪
  * zoomIn: 等比放大 
  * repeat: 平铺循环 
  * repeatX: 水平循环 
  * repeatY: 垂直循环 


### clip(cropper,shape)
对图像进行形状裁剪
```javascript
import Cropper from 'cropperjs';
this.editorObj.useCropper(Cropper,'circle');
```
* cropper : ``Object`` 裁剪的cropper插件对象
* shape : ``String`` 裁剪形状接受 ``square`` 和 ``circle``


### setColor(color)
接受一个颜色值，改变当前色块的颜色
```javascript
let color = 'red';
editorObj.setColor(color);
```


### blade.set(Number)
设置刀片的索引，可返回editor对象，可以在点击样板前先设置当前刀片的索引之后再设置样板,默认值为0
```javascript
editorObj.blade.set(1).setTemp(sampleInfo);
```
* Number : ``Number`` 刀片下标


### blade.get()
返回当前刀片的索引
```javascript
editorObj.blade.get();
```
* return : ``Number`` 刀片下标


### getPreviewDom()
获取预览区的DOM，用于查看大图渲染
```javascript
editorObj.getPreviewDom();
```
* return : ``String`` html dom


### clearEditor()
清空当前刀片的设计内容
```javascript
editorObj.clearPreview();
```


### nextStep()
调用历史记录的下一步
```javascript
editorObj.nextStep();
```

### preStep()
调用历史记录的上一步
```javascript
editorObj.preStep();
```

