/* eslint-disable */
// 获取 LocalStornag
const getLocalStornage = (name) => {
    return window.localStorage.getItem(name);
};
// 判断是否支持webp格式的图片
const util = {
    isSupportWebp: () => {
        const isSupportWebp = getLocalStornage('isSupportWebp');
        if (isSupportWebp) {
            return true;
        }
        try {
            const isSupportWebp =
                document
                    .createElement('canvas')
                    .toDataURL('image/webp', 0.5)
                    .indexOf('data:image/webp') === 0;
            setLocalStornage('isSupportWebp', true);
            return isSupportWebp;
        } catch (err) {
            return false;
        }
    }
}

class CanvasDistort {
    /**
     * 设置Canvas, 绘制在其内的Image, 将做为变形的素材.
     * @param {Canvas} canvas Canvas 实例
     * @param {number} imgStartX 绘图起始点X坐标
     * @param {number} imgStartY 绘图起始点Y坐标
     */
    constructor(
        canvasConfig = {width: 800, height: 800},
        imgStartX = 0,
        imgStartY = 0
    ) {
        const {width, height} = canvasConfig;
        // 创建canvas元素
        this.canvas = document.createElement('canvas');
        // this.canvas = new OffscreenCanvas(width, height);
        // 设置canvas的宽度和高度

        this.canvas.width = width;
        this.canvas.height = height;

        this.ctx = this.canvas.getContext('2d');
        this.baseX = 0; //基准线横坐标 旨在通过基准线计算偏移量
        this.baseY = 0;
        this.imgStartX = imgStartX;
        this.imgStartY = imgStartY;
        this.imgWidth = 0;
        this.imgHeight = 0;
        this.imgChangeObj = {};
        this.imgData = '';
        this.imgDataArr = {}; //存放图像原始数据与暂时数据
        this.imgAnimateDataObj = {}; //存储图像合成后各个图像imgData数据
    }

    /**
     * 弯曲Canvas中的图处
     * @param {string}} imgSrc 图片的src
     * @param {[]} ctrlNodes Bezier曲线的控制点坐标列表
     * @param {向哪个方向弯曲} type 弯曲方向, 可选择的值: 'col', 'row', 分别表示纵向弯曲, 横向弯曲
     * @param {(imgData)=>{})} callback 回调方法, 参数是弯曲之后图片的数据.
     *
     * 关于ctrlNodes的坐标设置:
     * 1. 当图片加载到canvas之后, 将设置X/Y方向的基线, 其中 X = imgStartX + imgWidth / 2, Y = imgStartY + imgHeight / 2
     * 2. ctrlNodes坐标x, y设置的范围:
     *     type = 'col', 满足 imgStartX <= x <= imgStartX + imgWidth, Y <= y <= Y + imgHeight / 2
     *     type = 'row', 满足 X <= x <= X + imgWidth / 2, imgStartY <= y <= imgStartY + imgHeight
     */
    distort = (imgSrc, ctrlNodes, type, callback) => {
        this.img = new Image();
        this.img.src = imgSrc;
        this.img.onload = () => {
            this.ctx.clearRect(
                this.imgStartX,
                this.imgStartY,
                this.imgWidth,
                this.imgHeight
            );
            this.imgWidth = this.img.width;
            this.imgHeight = this.img.height;
            this.ctx.drawImage(this.img, this.imgStartX, this.imgStartY);
            this.baseX = this.imgStartX + this.imgWidth / 2;
            this.baseY = this.imgStartY + this.imgHeight / 2;
            this.imgData = this.ctx.getImageData(
                this.imgStartX,
                this.imgStartY,
                this.imgWidth,
                this.imgHeight
            );
            this.imgDataArr.origin = [];
            this.imgData.data.forEach((item, index) => {
                this.imgDataArr.origin.push(item);
            });

            const imagedata = this.distortImg(ctrlNodes, type);
            callback(imagedata);
        };
    };

    distortImg = (ctrlNodes, type) => {
        this.imgAnimateDataObj = this.sliceImgData(ctrlNodes, type); //控制点转换为生成的扭曲图像重新赋值
        const canvasBg = document.createElement('canvas'); //离屏canvas 导出扭曲图片
        canvasBg.width = this.imgWidth;
        canvasBg.height = this.imgHeight;
        const bgctx = canvasBg.getContext('2d');
        let imgData = bgctx.getImageData(0, 0, this.imgWidth, this.imgHeight);
        this.imgAnimateDataObj.forEach((item, index) => {
            imgData.data[index] = item;
        });
        bgctx.putImageData(imgData, 0, 0);

        const isSupportWebp = util.isSupportWebp();
        const imgType = isSupportWebp ? 'webp' : 'png';

        const imgBse64 = canvasBg.toDataURL('image/' + imgType, 1);
        return imgBse64;
    };

    sliceImgData = (ctrlNodes, type) => {
        //对图像进行切分，分为row与col两种切分方式
        let arr = [];
        let bezierArr = [];
        let imgDataSlice = [];
        for (let i = 0; i < 1; i += 0.001) {
            bezierArr.push(this.bezier(ctrlNodes, i));
        }

        this.imgDataArr.origin.forEach((item, index) => {
            //移位后像素状态会产生变化需要重置
            this.imgData.data[index] = item;
        });

        bezierArr.forEach((obj, index) => {
            if (
                this.imgStartY < obj.y &&
                this.imgStartY + this.imgHeight > obj.y &&
                type === 'row'
            ) {
                let diffX = parseInt(obj.x - this.baseX, 10); //计算偏移量
                let dissY = parseInt(obj.y - this.imgStartY, 10);
                let rowNum = dissY;
                imgDataSlice = this.imgData.data.slice(
                    rowNum * this.imgWidth * 4,
                    rowNum * this.imgWidth * 4 + this.imgWidth * 4
                ); //按层切片
                for (var i = 0; i < Math.abs(diffX * 4); i++) {
                    imgDataSlice = this.arraymove(diffX, imgDataSlice);
                }
                this.imgChangeObj[rowNum] = imgDataSlice;
            } else if (
                this.baseX + this.imgWidth / 2 > obj.x &&
                this.baseX - this.imgWidth / 2 < obj.x &&
                type === 'col'
            ) {
                let diffX = parseInt(obj.x - (this.baseX - this.imgWidth / 2), 10); //计算偏移量
                let diffY = parseInt(obj.y - this.baseY, 10);
                let rowNum = diffX;
                this.imgChangeObj[rowNum] = {
                    diffX: diffX,
                    diffY: diffY,
                };
            }
        });

        if (type === 'row') {
            Object.keys(this.imgChangeObj).forEach((item, index) => {
                arr = arr.concat(Array.from(this.imgChangeObj[item]));
            });
        } else if (type === 'col') {
            for (let i = 0; i < this.imgWidth; i++) {
                imgDataSlice = [];
                for (let j = 0; j < this.imgHeight; j++) {
                    let index = j * this.imgWidth * 4 + i * 4;
                    let sliceArr = this.imgData.data.slice(index, index + 4);
                    imgDataSlice = imgDataSlice.concat(Array.from(sliceArr));
                }
                if (this.imgChangeObj[i]) {
                    for (let k = 0; k < Math.abs(this.imgChangeObj[i].diffY * 4); k++) {
                        imgDataSlice = this.arraymove(
                            this.imgChangeObj[i].diffY,
                            imgDataSlice
                        );
                    }
                    for (let p = 0; p < imgDataSlice.length / 4; p++) {
                        arr[p * this.imgWidth * 4 + i * 4] = imgDataSlice[p * 4];
                        arr[p * this.imgWidth * 4 + i * 4 + 1] = imgDataSlice[p * 4 + 1];
                        arr[p * this.imgWidth * 4 + i * 4 + 2] = imgDataSlice[p * 4 + 2];
                        arr[p * this.imgWidth * 4 + i * 4 + 3] = imgDataSlice[p * 4 + 3];
                    }
                }
            }
        }
        return arr;
    };

    arraymove = (type, arr) => {
        //切片后的数据进行移位
        let newArray = [];
        if (type > 0) {
            //右移
            let lastOne = arr[arr.length - 1];
            for (let i = 0; i < arr.length - 1; i++) {
                newArray[i + 1] = arr[i];
            }
            newArray[0] = lastOne;
        } else {
            let firstOne = arr[0];
            for (let i = 1; i < arr.length; i++) {
                newArray[i - 1] = arr[i];
            }
            newArray[arr.length - 1] = firstOne;
        }
        return newArray;
    };

    factorial = (num) => {
        //计算阶乘
        if (num <= 1) {
            return 1;
        } else {
            return num * this.factorial(num - 1);
        }
    };

    bezier = (arr, t) => {
        //通过各控制点与占比t计算当前贝塞尔曲线上的点坐标，坐标点不同是由t的值决定的。控制点数组均不变
        let x = 0,
            y = 0,
            n = arr.length - 1;

        arr.forEach((item, index) => {
            if (!index) {
                x += item.x * Math.pow(1 - t, n - index) * Math.pow(t, index);
                y += item.y * Math.pow(1 - t, n - index) * Math.pow(t, index);
            } else {
                x +=
                    (this.factorial(n) /
                        this.factorial(index) /
                        this.factorial(n - index)) *
                    item.x *
                    Math.pow(1 - t, n - index) *
                    Math.pow(t, index);
                y +=
                    (this.factorial(n) /
                        this.factorial(index) /
                        this.factorial(n - index)) *
                    item.y *
                    Math.pow(1 - t, n - index) *
                    Math.pow(t, index);
            }
        });
        return {x, y};
    };
}

// 2D 点位
class Point2D {
    constructor(x, y, u, v) {
        this.x = x;
        this.y = y;
        this.u = u;
        this.v = v;
    }

    clone() {
        return new Point2D(this.x, this.y, this.u, this.v);
    }
}

const DRAW_IMAGE_EXTEND_EX = 3;

// 2d 画图
class Vert2D {
    constructor(p0, p1, p2) {
        this.p0 = p0;
        this.p1 = p1;
        this.p2 = p2;
    }

    clone() {
        return new Vert2D(this.p0, this.p1, this.p2);
    }

    drawMeshLineToContext(plist, ctx) {
        let p0 = plist[this.p0],
            p1 = plist[this.p1],
            p2 = plist[this.p2];
        ctx.moveTo(p0.x, p0.y);
        ctx.lineTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.lineTo(p0.x, p0.y);
    }

    drawImageToContext(plist, img, ctx) {
        let p0 = plist[this.p0],
            p1 = plist[this.p1],
            p2 = plist[this.p2];
        Vert2D.drawImageToContextWithPoints(
            img,
            ctx,
            p0.x,
            p0.y,
            p1.x,
            p1.y,
            p2.x,
            p2.y,
            p0.u,
            p0.v,
            p1.u,
            p1.v,
            p2.u,
            p2.v
        );
    }

    static extendVert(x0, y0, x1, y1, x2, y2) {
        let x = 2 * x0 - x1 - x2,
            y = 2 * y0 - y1 - y2;
        let d = Math.sqrt(DRAW_IMAGE_EXTEND_EX / (x * x + y * y));
        return [x0 + x * d, y0 + y * d];
    }

    static drawImageToContextWithPoints(
        img,
        ctx,
        x0,
        y0,
        x1,
        y1,
        x2,
        y2,
        u0,
        v0,
        u1,
        v1,
        u2,
        v2
    ) {
        u0 *= img.width;
        u1 *= img.width;
        u2 *= img.width;
        v0 *= img.height;
        v1 *= img.height;
        v2 *= img.height;

        //fix gap in images
        let s0 = Vert2D.extendVert(x0, y0, x1, y1, x2, y2);
        let s1 = Vert2D.extendVert(x1, y1, x0, y0, x2, y2);
        let s2 = Vert2D.extendVert(x2, y2, x1, y1, x0, y0);
        //fix end

        ctx.beginPath();
        ctx.moveTo(s0[0], s0[1]);
        ctx.lineTo(s1[0], s1[1]);
        ctx.lineTo(s2[0], s2[1]);
        ctx.closePath();

        x1 -= x0;
        y1 -= y0;
        x2 -= x0;
        y2 -= y0;

        u1 -= u0;
        v1 -= v0;
        u2 -= u0;
        v2 -= v0;

        let det = 1 / (u1 * v2 - u2 * v1),
            a = (v2 * x1 - v1 * x2) * det,
            b = (v2 * y1 - v1 * y2) * det,
            c = (u1 * x2 - u2 * x1) * det,
            d = (u1 * y2 - u2 * y1) * det,
            e = x0 - a * u0 - c * v0,
            f = y0 - b * u0 - d * v0;

        ctx.save();
        ctx.transform(a, b, c, d, e, f);
        ctx.clip();
        ctx.drawImage(img, 0, 0);
        ctx.restore();
    }
}

class Mesh2D {
    constructor() {
        this.points = [];
        this.verts = [];
    }

    clone() {
        let n = new Mesh2D();
        for (let i = 0; i < this.points.length; i++) {
            n.points[i] = this.points[i].clone();
        }
        for (let i = 0; i < this.verts.length; i++) {
            n.verts[i] = this.verts[i];
        }
        return n;
    }

    move(x, y) {
        for (let i = 0; i < this.points.length; i++) {
            this.points[i].x += x;
            this.points[i].y += y;
        }
    }

    changeByMatrix4(te) {
        for (let i = 0; i < this.points.length; i++) {
            this.points[i].changeByMatrix4(te);
        }
    }

    drawImageToContext(img, ctx) {
        for (let i = 0; i < this.verts.length; i++) {
            this.verts[i].drawImageToContext(this.points, img, ctx);
        }
    }

    static createMapMesh(width, height, divW, divH) {
        let m = new Mesh2D();
        let widthSingle = width / divW,
            heightSingle = height / divH;
        let uSingle = 1 / divW,
            vSingel = 1 / divH;
        for (let i = 0; i <= divH; i++) {
            for (let j = 0; j <= divW; j++) {
                m.points.push(
                    new Point2D(
                        j * widthSingle,
                        i * heightSingle,
                        j * uSingle,
                        i * vSingel
                    )
                );
            }
        }
        for (let i = 0; i < divH; i++) {
            for (let j = 0; j < divW; j++) {
                let startPoint = (divW + 1) * i + j;
                m.verts.push(
                    new Vert2D(startPoint + 1, startPoint, startPoint + divW + 1)
                );
                m.verts.push(
                    new Vert2D(
                        startPoint + divW + 1,
                        startPoint + divW + 2,
                        startPoint + 1
                    )
                );
            }
        }
        return m;
    }
}

// 曲线算法
class Mumeric {
    static dim(x) {
        let y, z;
        if (typeof x === 'object') {
            y = x[0];
            if (typeof y === 'object') {
                z = y[0];
                if (typeof z === 'object') {
                    return numeric._dim(x);
                }
                return [x.length, y.length];
            }
            return [x.length];
        }
        return [];
    }

    static _foreach2(x, s, k, f) {
        if (k === s.length - 1) {
            return f(x);
        }
        let i,
            n = s[k],
            ret = Array(n);
        for (i = n - 1; i >= 0; i--) {
            ret[i] = this._foreach2(x[i], s, k + 1, f);
        }
        return ret;
    }

    static cloneV(x) {
        let _n = x.length;
        let i,
            ret = Array(_n);

        for (i = _n - 1; i !== -1; --i) {
            ret[i] = x[i];
        }
        return ret;
    }

    static clone(x) {
        if (typeof x !== 'object') return x;
        let V = this.cloneV;
        let s = this.dim(x);
        return this._foreach2(x, s, 0, V);
    }

    static diag(d) {
        let i,
            i1,
            j,
            n = d.length,
            A = Array(n),
            Ai;
        for (i = n - 1; i >= 0; i--) {
            Ai = Array(n);
            i1 = i + 2;
            for (j = n - 1; j >= i1; j -= 2) {
                Ai[j] = 0;
                Ai[j - 1] = 0;
            }
            if (j > i) {
                Ai[j] = 0;
            }
            Ai[i] = d[i];
            for (j = i - 1; j >= 1; j -= 2) {
                Ai[j] = 0;
                Ai[j - 1] = 0;
            }
            if (j === 0) {
                Ai[0] = 0;
            }
            A[i] = Ai;
        }
        return A;
    }

    static rep(s, v, k) {
        if (typeof k === 'undefined') {
            k = 0;
        }
        let n = s[k],
            ret = Array(n),
            i;
        if (k === s.length - 1) {
            for (i = n - 2; i >= 0; i -= 2) {
                ret[i + 1] = v;
                ret[i] = v;
            }
            if (i === -1) {
                ret[0] = v;
            }
            return ret;
        }
        for (i = n - 1; i >= 0; i--) {
            ret[i] = this.rep(s, v, k + 1);
        }
        return ret;
    }

    static identity(n) {
        return this.diag(this.rep([n], 1));
    }

    static inv(a) {
        let s = this.dim(a),
            abs = Math.abs,
            m = s[0],
            n = s[1];
        let A = this.clone(a),
            Ai,
            Aj;
        let I = this.identity(m),
            Ii,
            Ij;
        let i, j, k, x;
        for (j = 0; j < n; ++j) {
            let i0 = -1;
            let v0 = -1;
            for (i = j; i !== m; ++i) {
                k = abs(A[i][j]);
                if (k > v0) {
                    i0 = i;
                    v0 = k;
                }
            }
            Aj = A[i0];
            A[i0] = A[j];
            A[j] = Aj;
            Ij = I[i0];
            I[i0] = I[j];
            I[j] = Ij;
            x = Aj[j];
            for (k = j; k !== n; ++k) Aj[k] /= x;
            for (k = n - 1; k !== -1; --k) Ij[k] /= x;
            for (i = m - 1; i !== -1; --i) {
                if (i !== j) {
                    Ai = A[i];
                    Ii = I[i];
                    x = Ai[j];
                    for (k = j + 1; k !== n; ++k) Ai[k] -= Aj[k] * x;
                    for (k = n - 1; k > 0; --k) {
                        Ii[k] -= Ij[k] * x;
                        --k;
                        Ii[k] -= Ij[k] * x;
                    }
                    if (k === 0) Ii[0] -= Ij[0] * x;
                }
            }
        }
        return I;
    }

    static dotMMsmall(x, y) {
        let i, j, k, p, q, r, ret, foo, bar, woo, i0;
        p = x.length;
        q = y.length;
        r = y[0].length;
        ret = Array(p);
        for (i = p - 1; i >= 0; i--) {
            foo = Array(r);
            bar = x[i];
            for (k = r - 1; k >= 0; k--) {
                woo = bar[q - 1] * y[q - 1][k];
                for (j = q - 2; j >= 1; j -= 2) {
                    i0 = j - 1;
                    woo += bar[j] * y[j][k] + bar[i0] * y[i0][k];
                }
                if (j === 0) {
                    woo += bar[0] * y[0][k];
                }
                foo[k] = woo;
            }
            ret[i] = foo;
        }
        return ret;
    }

    static dotMV(x, y) {
        let p = x.length,
            i;
        let ret = Array(p);
        // dotVV = this.dotVV;
        for (i = p - 1; i >= 0; i--) {
            ret[i] = this.dotVV(x[i], y);
        }
        return ret;
    }

    static dotVV(x, y) {
        let i,
            n = x.length,
            i1,
            ret = x[n - 1] * y[n - 1];
        for (i = n - 2; i >= 1; i -= 2) {
            i1 = i - 1;
            ret += x[i] * y[i] + x[i1] * y[i1];
        }
        if (i === 0) {
            ret += x[0] * y[0];
        }
        return ret;
    }

    static transpose(x) {
        let i,
            j,
            m = x.length,
            n = x[0].length,
            ret = Array(n),
            A0,
            A1,
            Bj;
        for (j = 0; j < n; j++) ret[j] = Array(m);
        for (i = m - 1; i >= 1; i -= 2) {
            A1 = x[i];
            A0 = x[i - 1];
            for (j = n - 1; j >= 1; --j) {
                Bj = ret[j];
                Bj[i] = A1[j];
                Bj[i - 1] = A0[j];
                --j;
                Bj = ret[j];
                Bj[i] = A1[j];
                Bj[i - 1] = A0[j];
            }
            if (j === 0) {
                Bj = ret[0];
                Bj[i] = A1[0];
                Bj[i - 1] = A0[0];
            }
        }
        if (i === 0) {
            A0 = x[0];
            for (j = n - 1; j >= 1; --j) {
                ret[j][0] = A0[j];
                --j;
                ret[j][0] = A0[j];
            }
            if (j === 0) {
                ret[0][0] = A0[0];
            }
        }
        return ret;
    }
}

// 透视法转换图像
class Perspective {
    constructor(srcPts, dstPts) {
        this.srcPts = srcPts;
        this.dstPts = dstPts;
        this.coeffs = this.getNormalizationCoefficients(
            this.srcPts,
            this.dstPts,
            false
        );
        this.coeffsInv = this.getNormalizationCoefficients(
            this.srcPts,
            this.dstPts,
            true
        );
    }

    getNormalizationCoefficients(srcPts, dstPts, isInverse) {
        if (isInverse) {
            let tmp = dstPts;
            dstPts = srcPts;
            srcPts = tmp;
        }
        let r1 = [
            srcPts[0],
            srcPts[1],
            1,
            0,
            0,
            0,
            -1 * dstPts[0] * srcPts[0],
            -1 * dstPts[0] * srcPts[1],
        ];
        let r2 = [
            0,
            0,
            0,
            srcPts[0],
            srcPts[1],
            1,
            -1 * dstPts[1] * srcPts[0],
            -1 * dstPts[1] * srcPts[1],
        ];
        let r3 = [
            srcPts[2],
            srcPts[3],
            1,
            0,
            0,
            0,
            -1 * dstPts[2] * srcPts[2],
            -1 * dstPts[2] * srcPts[3],
        ];
        let r4 = [
            0,
            0,
            0,
            srcPts[2],
            srcPts[3],
            1,
            -1 * dstPts[3] * srcPts[2],
            -1 * dstPts[3] * srcPts[3],
        ];
        let r5 = [
            srcPts[4],
            srcPts[5],
            1,
            0,
            0,
            0,
            -1 * dstPts[4] * srcPts[4],
            -1 * dstPts[4] * srcPts[5],
        ];
        let r6 = [
            0,
            0,
            0,
            srcPts[4],
            srcPts[5],
            1,
            -1 * dstPts[5] * srcPts[4],
            -1 * dstPts[5] * srcPts[5],
        ];
        let r7 = [
            srcPts[6],
            srcPts[7],
            1,
            0,
            0,
            0,
            -1 * dstPts[6] * srcPts[6],
            -1 * dstPts[6] * srcPts[7],
        ];
        let r8 = [
            0,
            0,
            0,
            srcPts[6],
            srcPts[7],
            1,
            -1 * dstPts[7] * srcPts[6],
            -1 * dstPts[7] * srcPts[7],
        ];

        let matA = [r1, r2, r3, r4, r5, r6, r7, r8];
        let matB = dstPts;
        let matC;

        try {
            matC = Mumeric.inv(Mumeric.dotMMsmall(Mumeric.transpose(matA), matA));
        } catch (e) {
            return [1, 0, 0, 0, 1, 0, 0, 0];
        }

        let matD = Mumeric.dotMMsmall(matC, Mumeric.transpose(matA));
        let matX = Mumeric.dotMV(matD, matB);
        for (let i = 0; i < matX.length; i++) {
            matX[i] = Math.round(matX[i] * 10000000000) / 10000000000;
        }
        matX[8] = 1;

        return matX;
    }

    static transform(x, y) {
        let coordinates = [];
        coordinates[0] =
            (this.coeffs[0] * x + this.coeffs[1] * y + this.coeffs[2]) /
            (this.coeffs[6] * x + this.coeffs[7] * y + 1);
        coordinates[1] =
            (this.coeffs[3] * x + this.coeffs[4] * y + this.coeffs[5]) /
            (this.coeffs[6] * x + this.coeffs[7] * y + 1);
        return coordinates;
    }

    static transformFun(x, y, coeffs) {
        let coordinates = [];
        coordinates[0] =
            (coeffs[0] * x + coeffs[1] * y + coeffs[2]) /
            (coeffs[6] * x + coeffs[7] * y + 1);
        coordinates[1] =
            (coeffs[3] * x + coeffs[4] * y + coeffs[5]) /
            (coeffs[6] * x + coeffs[7] * y + 1);
        return coordinates;
    }

    static transformInverse(x, y) {
        let coordinates = [];
        coordinates[0] =
            (this.coeffsInv[0] * x + this.coeffsInv[1] * y + this.coeffsInv[2]) /
            (this.coeffsInv[6] * x + this.coeffsInv[7] * y + 1);
        coordinates[1] =
            (this.coeffsInv[3] * x + this.coeffsInv[4] * y + this.coeffsInv[5]) /
            (this.coeffsInv[6] * x + this.coeffsInv[7] * y + 1);
        return coordinates;
    }
}

// canvas画布管理器
class CanvasManager {
    constructor(canvasConfig = {width: 600, height: 600}) {
        const {width, height} = canvasConfig;
        // 创建canvas元素
        this.canvas = document.createElement('canvas');
        // console.log(this.canvas, 'canvas')
        // 设置canvas的宽度和高度

        this.canvas.width = width;
        this.canvas.height = height;

        // 设置载入图片时的宽度和高度,
        // 一般设置为铺满整个canvas,
        // 如果不铺满, 产生正视图时会变形
        const itemWidth = width;
        const itemHeight = height;

        // 设置ViewPort的相对于Canvas的偏移
        const viewPortOffsetX = 0;
        const viewPortOffsetY = 0.5;

        // 设置ViewPort的宽,高
        // 一般设置为正向面对canvas, 且贴近于canvas.
        // 如果设置的宽/高大于canvas的宽高, 生成的图像会缩小
        // 反之, 图像会放大并产生剪切效果
        const viewPortWidth = width;
        const viewPortHeight = height;

        const viewPortCoordinates = [
            {
                x: viewPortOffsetX,
                y: viewPortOffsetY,
            },
            {
                x: viewPortOffsetX + viewPortWidth,
                y: viewPortOffsetY,
            },
            {
                x: viewPortOffsetX + viewPortWidth,
                y: viewPortOffsetY + viewPortHeight,
            },
            {
                x: viewPortOffsetX,
                y: viewPortOffsetY + viewPortHeight,
            },
        ];

        this.viewPort = viewPortCoordinates;
        this.events = {};
        this.ctx = this.canvas.getContext('2d');
        // 填充背景为透明
        this.ctx.fillStyle = 'rgba(255, 255, 255, 0.0)';
        // 绘制矩形区域
        this.ctx.fillRect(0, 0, width, height);

        const centerX = width / 2;
        const centerY = height / 2;
        this.ctx.translate(centerX, centerY);
        this.setImageDimension(width, width);
    }

    /**
     * 设置图片的宽/高信息
     * @param {number} imgW 原始图片宽度
     * @param {number} imgH 原始图片高度
     * @param {number} divW 宽度分成几份, 默认16, 分数越高, 对系统绘图性能越高
     * @param {number} divH 高度分成几份, 默认16
     */
    setImageDimension = (imgW, imgH, divW = 10, divH = 10) => {
        this.model = Mesh2D.createMapMesh(imgW, imgH, divW, divH);
        return this;
    };
    /**
     * 设置图片的[左上, 右上, 右下, 左下]四个角的坐标
     * @param {[Object]} coordinates 图片四个角的坐标[{x: -100, y: -100}, {x: 100, y: -100}, {x: 100, y: 100}, {x: -100, y: 100}]
     */
    setImageCoordinates = () => {
        this.srcPoints = this.viewPort;
        return this;
    };
    /**
     * 设置输出区域的四个角的坐标
     * @param {[Object]} coordinates 输出区域四角的坐标
     */
    setOutputCoordinates = (coordinates) => {
        this.dstPoints = coordinates;
        return this;
    };
    /**
     * 移动图像, 在已经设置好其它参数, 并已经draw之后调用.
     * @param {number} offsetX x轴偏移量
     * @param {number} offsetY y轴偏移量
     */
    translate = (offsetX, offsetY) => {
        this.model.move(offsetX, offsetY);
        return this;
    };
    /**
     * 生成透视图
     */
    perspective = () => {
        let n = this.srcPoints;
        let l = this.dstPoints;
        // console.log(this.srcPoints, this.dstPoints,'==================')
        let e = [n[0].x, n[0].y, n[1].x, n[1].y, n[2].x, n[2].y, n[3].x, n[3].y];
        let t = [l[0].x, l[0].y, l[1].x, l[1].y, l[2].x, l[2].y, l[3].x, l[3].y];
        // console.log(e,t,'----------------------------')
        let perspT = new Perspective(e, t);

        const {coeffs} = perspT;

        for (let e = 0; e < this.model.points.length; e++) {
            let t = Perspective.transformFun(
                this.model.points[e].x,
                this.model.points[e].y,
                coeffs
            );

            (this.model.points[e].x = t[0]), (this.model.points[e].y = t[1]);
        }

        return this;
        return this;
    };
    /**
     * 清空绘制内容
     */
    clear = () => {
        this.ctx.save();
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.ctx.restore();
        return this;
    };
    genImg = (canvasArea, imgBase64, moveX, moveY) => {
        return new Promise((resolve) => {
            this.setImageCoordinates();
            this.setOutputCoordinates(canvasArea);
            this.perspective();
            this.clear();

            const img = new Image();
            if (imgBase64.search('base64') != -1) {
                img.src = imgBase64;
            } else if (imgBase64.search('x-oss-process') != -1) {
                img.src = imgBase64.split('?')[0] + '?time=' + new Date().getTime();
            } else {
                img.src = imgBase64 + '?time=' + new Date().getTime();
            }
            img.setAttribute('crossOrigin', 'Anonymous');

            const isSupportWebp = util.isSupportWebp();
            const type = isSupportWebp ? 'webp' : 'png';

            img.onload = () => {
                // console.log(left,top, 'this.model.points');
                this.model.drawImageToContext(img, this.ctx);
                const base64 = this.canvas.toDataURL('image/' + type, 1);
                resolve(base64);
            };
        });
    };
}

// 坐标点转换成贴图的区域
const stringToArray = (arr, canvasOptions) => {
    const arrList = [];
    let {width, height} = canvasOptions;
    const defaultWidth = 800;
    const defaultHeight = 800;
    let scaleX = 1;
    let scaleY = 1;
    if (width) {
        scaleX = width / defaultWidth;
        scaleY = height / defaultHeight;
    }

    for (let i = 0; i < arr.length; i++) {
        const item = {
            x: arr[i].split(',')[0] * scaleX,
            y: arr[i].split(',')[1] * scaleY,
        };
        arrList.push(item);
    }
    return arrList;
};

// 历史记录管理
class HistoryManager {
    constructor() {
        this.HISTORY = {
            color:{

            },
            pattern:{

            },
            text:{

            }
        }
    }
    find(obj){
        let type = this.getObjType(obj); // 获取当前对象类型
        let id = obj.getAttribute(type + 'id'); // id是由当前对象类型加id拼接的。
        if(!this.HISTORY[type][type+id]) return [];
        return this.HISTORY[type][type+id][this.HISTORY[type][type+id].length - 1];
    }
    delete(obj){

    }
    add(obj){
        let type = this.getObjType(obj);
        let info = this.getInfo(obj, type);
        // 先判断记录里有没有和这个对象的记录,
        let key = type+info.id;
        console.log(type,info,key)
        if(!this.HISTORY[type][key]){
            Object.assign(this.HISTORY[type],{[key]: []});
        }

        // 如果超过10条则需要删除第一条
        if(this.HISTORY[type][key] && this.HISTORY[type][key].length === 10){
            this.HISTORY[type][key].shift();
        }
        this.HISTORY[type][key].push(info);
        console.log(this.HISTORY[type][key])

    }
    // 获取要存储的信息
    getInfo(obj, type){
        const info = {
            type: type,
            id: obj.getAttribute('colorid') || obj.getAttribute('patternid'),
            width: obj.clientWidth,
            height: obj.clientHeight,
            left: obj.offsetLeft,
            top: obj.offsetTop,
            transform: obj.style.transform,
            zIndex: obj.style.zIndex,
        };
        if(type === 'pattern'){
            info.src = obj.src;
        }
        if(type === 'color'){
            info.colorType = obj.getAttribute('colortype');
        }
        return info;
    }
    // 获取当前对象的类型
    getObjType(obj){
        return obj.getAttribute('data-type')
    }
}
const point2D = Point2D;
const vert2D = Vert2D;
const mesh2D = Mesh2D;
const perspective = Perspective;
const canvasManager = CanvasManager;
const discord = CanvasDistort;
/**
 * IYWDesigner
 * @constructor
 * @param {String} id - 设计器容器ID.
 * @param {String} previewBox - 预览区ID.
 * @param {Number} exportSize - 效果图的尺寸.
 * @param {Number} maxItemNum - 设计元素的最大数量.
 * @author Jerry date 2022/05/03
 */
export default class IYWDesigner {
    constructor({id, previewBox, selectedInfo = {info: {}}, exportSize = 800, maxItemNum = 6}) {
        this.selectedInfo = selectedInfo; // 当前操作对象的详情会通过这个对象返回到调用者
        this.exportSize = exportSize; // 导出效果图的尺寸
        this.exportObj = null; // 导出效果图的盒子
        this.canvasSize = 600;
        this.currentBlade = 0;
        this.frame = null; // 操作框
        this.id = id; // 操作区域的父级元素
        this.previewBox = previewBox; // 预览区的盒子
        this.image = null // 当前选中的图片
        this.editorBox = null; // 生成编辑区域dom
        this.disX = null; // 鼠标移动的坐标
        this.disY = null;
        this.currentAction = ''; //
        this.cropper = null; //裁剪插件
        this.previewImgContainer = null;
        this.onMoveObj = null; // 当前操作的对象
        this.defaultSample = null; // 样版信息
        this.isDrag = false; // 控制图片拖拽
        this.isHold = false; // 控制图片旋转
        this.isTopHold = false; // 控制编辑器操作点上
        this.isBottomHold = false;
        this.isLeftpHold = false;
        this.isRightpHold = false;
        this.historyManage = new HistoryManager();
        this.angel = 0; // 旋转的角度
        this.objDefParams = {};
        this.defaultLeft = 0; // 编辑器框的位置
        this.defaultTop = 0;
        this.controlElIndex = null; // 矩阵用
        this.elementsParent = null; // 设计器内部所有可操作元素的父类
        this.layers = [];
        this.maxItemNum = maxItemNum; // 设计区内元素的最大长度
        this.pointA = {
            x: 0,
            y: 0,
        };
        this.pointB = {
            x: 0,
            y: 0,
        };
        this.pointC = {
            x: 0,
            y: 0
        };
        this.allAngel = 0; // 鼠标角度
        this.count = 0;
        this.colorConfig = {
            size: 300,
            data: []
        }
        this.frameParentName = "editor-ctrl"; // 编辑框操作节点组的父节点
        this.designContents = { // 当前设计内容
            bladeList: [
                {
                    bladeImg: "https://chdesign.oss-cn-shanghai.aliyuncs.com/test/CUMS/OA/20211111/163662509799073626.png",
                    bladeName: "刀片1",
                    contents: {
                        images: null,
                        colorPiece: [],
                        text: []
                    }
                }
            ],
            template: null
        }
        // 刀片操作方法集合
        this.blade = {
            set: (val) => {
                this.currentBlade = val;
                return this;
            },
            get: () => {
                return this.currentBlade;
            },
            replace: () => {

            },
            find: () => {

            },
            getCurrentBlade: (badeIndex, colorIndex) => {
                return this.defaultSample.baseImgInfoList[this.currentBlade].bladeInfoList[colorIndex]
            }
        }
        // 要生成的操作节点
        this.framePoints = [
            "top", "top-left", "top-right", "left", "left-bottom", "bottom", "bottom-right", "right"
        ]
    }

    // 获取预览区的节点
    getPreviewDom(){
        return this.previewImgContainer.innerHTML;
    }

    showFrame() {
        this.frame.style.opacity = '1';
    }

    hideFrame() {
        this.frame.style.opacity = '0';
    }

    // 设置刀片下标
    setCurrentBlade(index) {
        this.currentBlade = index;
    }

    // 生成单个合成div
    createExportSingle() {
        let canvas = this.createElement({type: 'div', id: 'compositeCan'});
        canvas.style.width = this.exportSize + 'px';
        canvas.style.height = this.exportSize + 'px';
        this.insertElement(document.getElementsByTagName('body')[0], canvas);
    }

    // 使用裁剪插件
    useCropper(copper, shape = 'square') {
        this.createDialogForCropper(shape);
        const img = document.getElementById('imgForClip');
        if (copper) {
            this.cropper = new copper(img, {
                aspectRatio: 1,
                viewMode: 1,
                ready: function () {
                    if (shape === 'circle') {
                        document.getElementsByClassName('cropper-crop-box')[0].style.borderRadius = '50%';
                        document.getElementsByClassName('cropper-crop-box')[0].style.overflow = 'hidden';
                        document.getElementsByClassName('cropper-face')[0].style.borderRadius = '50%';
                        document.getElementsByClassName('cropper-view-box')[0].style.borderRadius = '50%';
                    }
                }
            });
        }
    }

    // 裁剪成圆形
    getRoundedCanvas(sourceCanvas) {
        let canvas = document.createElement('canvas');
        let context = canvas.getContext('2d');
        let width = sourceCanvas.width;
        let height = sourceCanvas.height;
        canvas.width = width;
        canvas.height = height;
        context.imageSmoothingEnabled = true;
        context.drawImage(sourceCanvas, 0, 0, width, height);
        context.globalCompositeOperation = 'destination-in';
        context.beginPath();
        context.arc(width / 2, height / 2, Math.min(width, height) / 2, 0, 2 * Math.PI, true);
        context.fill();
        return canvas;
    }

    // 生成一个裁剪的弹窗
    createDialogForCropper(shape) {
        const dialog = this.createElement({type: 'div', className: 'iyw-dialog'});// 裁剪弹窗
        let diaLogBody = this.createElement({type: 'div', className: 'dialog-body'}); // 弹窗内容父级元素
        let diaLogImg = this.createElement({type: 'img', id: 'imgForClip'}); // 弹窗图片容器
        let clipContainer = this.createElement({type: 'div', className: 'clip-container'}) // 弹窗的裁剪容器
        let dialogFoot = this.createElement({type: 'div', className: 'clip-footer'})
        let confirmButton = this.createElement({type: 'button', className: 'clip-confirm'}); // 按钮
        let cancelButton = this.createElement({type: 'button', className: 'clip-cancel'}); // 按钮
        confirmButton.innerText = '确定';
        cancelButton.innerText = '取消';
        // 关闭裁剪窗口销毁cropper实例
        const closeDialog = () => {
            dialog.remove();
            this.cropper = null;
        }
        // 确认裁剪
        confirmButton.onclick = () => {
            const _index = this.currentBlade;
            let clippedImg = '';
            // 如果是圆形要调用裁剪圆形的方法进行裁剪
            if (shape === 'circle') {
                clippedImg = this.getRoundedCanvas(this.cropper.getCroppedCanvas()).toDataURL('image/png');
            } else {
                clippedImg = this.cropper.getCroppedCanvas().toDataURL('image/png')
            }

            // 进行处理的图片要添加到原图所属对象当中去。
            Object.assign(this.designContents.bladeList[_index].contents.images, {canvasImg: clippedImg});

            // 调用添加图片的方法
            this.addImage(this.designContents.bladeList[_index].contents.images);
            closeDialog();
        }
        // 取消裁剪
        cancelButton.onclick = () => {
            closeDialog();
        }
        diaLogImg.src = this.onMoveObj.src;
        clipContainer.appendChild(diaLogImg);
        diaLogBody.appendChild(clipContainer);
        dialogFoot.appendChild(cancelButton);
        dialogFoot.appendChild(confirmButton);
        diaLogBody.appendChild(dialogFoot);
        dialog.appendChild(diaLogBody);
        this.insertElement(document.getElementById('app'), dialog);
    }

    // 保存对象操作后的最后坐标和尺寸
    saveLatestInfo() {
        const _this = this;
        _this.objDefParams.left = _this.onMoveObj.offsetLeft; // 存储下来为自由变换对象宽高
        _this.objDefParams.top = _this.onMoveObj.offsetTop;
        _this.objDefParams.width = _this.onMoveObj.clientWidth; // 存储下来为自由变换对象宽高
        _this.objDefParams.height = _this.onMoveObj.clientHeight;
    }

    // 添加鼠标监听
    addEventsForEditor() {
        // 处理有色块的情况下， 删除所有色块最后的选择框的尺寸以及坐标还是最后一个色块的问题。
        const fixColorDelete = () => {
            // 如果聚焦的对象为null 先检查是不是色块都删除了
            if (!_this.colorConfig.data.length) {
                // 如果没有色块的，但是图片对象在而且有图片则把聚焦对象设置成图片
                if (_this.hasEditImage()) {
                    _this.onMoveObj = _this.hasEditImage();
                    let {offsetLeft, offsetTop} = _this.onMoveObj;
                    _this.asyncPositionForFrame(_this.onMoveObj, offsetLeft, offsetTop);
                }
            }
        }

        // 绑定拖拽事件
        let _this = this;
        const stopMove = (e) => {
            // 鼠标松开取消所有监听事件
            const mouseDownHandle =
                {
                    isDrag: 'onMove',
                    isHold: 'onRotate',
                    isTopHold: 'topOnMove',
                    isLeftHold: 'leftOnMove',
                    isRightHold: 'rightOnMove',
                    isBottomHold: 'bottomOnMove'
                };
            Object.keys(mouseDownHandle).forEach(item => {
                _this[item] = false;
            })
            if (_this && _this.onMoveObj) {
                _this.saveLatestInfo();
                _this.onMoveObj.style.cursor = "default";
                _this.frame.style.cursor = "default";
            }

            // 检查有没有旋转
            if(_this.currentAction === 'rotate'){
                _this.genImgForAction('rotate');
            }
            // this.frame.style.pointerEvents = 'none';
        }
        this.editorBox.onmousedown = function (e) {
            e.preventDefault();
            e.stopPropagation();
            if (!_this.elementsParent.childNodes.length) {
                _this.onMoveObj = null;
                fixColorDelete();
                _this.hideFrame();
                _this.getCanvasSize();
                return;
            }
            _this.isDrag = true;
            _this.frame.style.pointerEvents = 'auto';
            // _this.asyncPositionForFrame(_this.onMoveObj, _this.onMoveObj.offsetLeft,_this.onMoveObj.offsetTop);
            _this.showFrame();
            if (e.target.id !== 'elementsParent' && e.target.id !== _this.id && e.target.nodeName !== 'SPAN') {
                _this.focusSelectedItem(e.target); // 把当前聚焦的对象传过去保存起来
                if (!_this.onMoveObj) {
                    _this.hideFrame();
                    _this.isDrag = false;
                    return;
                }
                _this.isDrag = true;
                _this.frame.style.cursor = "move";
                _this.onMoveObj.style.cursor = 'move';
                e = e || window.event; // 兼容IE写法
                _this.isDrag = true;
                // 获取鼠标位置
                let x = e.clientX,
                    y = e.clientY;

                // 鼠标相对图片的位置
                _this.disX = x - _this.onMoveObj.offsetLeft;
                _this.disY = y - _this.onMoveObj.offsetTop;
                // _this.addMouseEvent();
                // console.log(e, '鼠标按下了')
            } else {
                _this.isDrag = false;
            }

        }
        document.onmousemove = function (e) {
            e.preventDefault();
            e.stopPropagation();
            // 第一次保存原始图片的坐标和尺寸
            if (!_this.objDefParams && _this.onMoveObj) {
            }
            // 事件集合，每个鼠标移动绑定的事件都有一个值来对应。
            const mouseDownHandle =
                {
                    isDrag: 'onMove',
                    isHold: 'onRotate',
                    isTopHold: 'topOnMove',
                    isLeftHold: 'leftOnMove',
                    isRightHold: 'rightOnMove',
                    isBottomHold: 'bottomOnMove'
                };
            Object.keys(mouseDownHandle).forEach(item => {
                if (_this[item]) {
                    let runMethod = mouseDownHandle[item];
                    _this[runMethod](e);
                }
            })
        };
        // 鼠标起来后解除图片移动的控制
        document.onmouseup = function (e) {
            if (_this.onMoveObj) {
                _this.saveLatestInfo();
                e.preventDefault();
                e.stopPropagation();
                _this.showAllColorCtl('block');
                stopMove(e);
            }
            _this.editorBox.onmouseup = function (e) {
                // 获取上一步的操作记录和现在的信息匹配，如果没有发生改变则不执行同步动作
                let lastStep = _this.historyManage.find(_this.onMoveObj);
                let type = _this.getElementType();
                let currentData = _this.historyManage.getInfo(_this.onMoveObj,type);
                _this.selectedInfo.info = _this.createObjInfo();
                if(JSON.stringify(lastStep) != JSON.stringify(currentData)){
                    _this.loadingPreview.style.display = 'block';
                    _this.asyncPosition(_this.previewBox);
                }
            }
        };
    }

    // 初始化编辑框
    init() {
        if (!this.id) {
            alert("必须传入一个有效的元素ID");
            return;
        }
        if (!document.getElementById(this.id)) {
            alert("ID 元素不存在");
            return;
        }
        if (this.previewBox) {
            // 保存预览区域的dom
            this.previewBox = document.getElementById(this.previewBox);
            let loadingBoxParent = this.createElement({type:'div',className: 'loading-box'});
            let loadingContext = this.createElement({type: 'div', className:'loader'});
            loadingBoxParent.appendChild(loadingContext);
            this.insertElement(this.previewBox, loadingBoxParent);
            this.loadingPreview = loadingBoxParent;
            loadingBoxParent.style.display = 'none';
        }
        // 保存编辑器的dom
        this.editorBox = document.getElementById(this.id);
        this.editorBox.style.pointerEvents = 'auto';

        const EDITOR = this.editorBox;
        this.getCanvasSize();
        // 判断有没有移动操作框
        if (!this.defaultLeft && !this.defaultTop) {

            // 创建外层的移动操作框
            const move = this.createElement({type: 'div', id: 'move'});
            this.insertElement(EDITOR, move);
            this.frame = move;

            // 创建设计器内部所有可操作元素的父类,所有添加的图片和色块都会放入到这个元素下面
            this.elementsParent = this.createElement({type: 'div', id: 'elementsParent'});
            this.insertElement(EDITOR, this.elementsParent);

            // 创建外框的8个点
            this.createFramePoints(this.frame);

            // 创建合成效果图的隐藏dom
            this.createExportBox();
        }
        // 操作框默认不显示

        this.createExportSingle();
        this.hideFrame();
        this.addEventsForEditor();

    }

    // 创建合成的DIV盒子
    createExportBox() {
        if (!document.getElementById('element-for-export')) {
            let elementForExport = this.createElement({id: 'element-for-export', type: 'div'});
            // 设置导出图的尺寸
            elementForExport.style.width = this.exportSize + 'px';
            elementForExport.style.height = this.exportSize + 'px';
            elementForExport.style.position = 'absolute';
            elementForExport.style.zIndex = '-90';
            elementForExport.style.right = '0';
            elementForExport.style.bottom = '0';
            elementForExport.style.perspective = '1200px'
            // elementForExport.style.display = 'none';

            elementForExport.style.pointerEvents = 'none';
            let body = document.getElementsByTagName('body')[0];
            this.insertElement(body, elementForExport);
            this.exportObj = elementForExport;
            this.insertElement(this.exportObj, this.createElement({type: 'div', className: 'preview-container'}));
        }
    }

    // 创建元素
    /**
     * @constructor
     * @param {String} type - 要创建的标签类型 如 span div.
     * @param {String} id - 元素的ID (非必填).
     * @param {String} className - 元素的class (非必填).
     * @param {String} content - 元素内的内容 (非必填).
     * @param {Function} event - 元素的事件 (非必填).
     */
    createElement({type = '', id = '', className = ''}) {
        let Element = document.createElement(type);
        if (id) Element.id = id;
        if (className) Element.className = className;
        return Element;
    }

    // 添加元素到父级盒子, 默认全部在传入的id的内部添加元素
    insertElement(parent, element) {
        if (!parent || !element) return;
        parent.appendChild(
            element
        )
    }

    // 设置样版
    setTemp(sampleInfo) {
        this.defaultSample = {...sampleInfo}; // 保存当前样版信息
        const bladeImg = this.findElementByClass(this.editorBox, 'bladeImg');
        if(!sampleInfo.bladeImgInfoList[this.currentBlade]){
            throw new Error('blade info is empty, found tId == '+ sampleInfo.tId)
        }
        let bladeImage = sampleInfo.bladeImgInfoList[this.currentBlade].bladeImg; // 刀片图
        // 设置刀片图，没有就生成一个img元素，有就直接替换src
        if (bladeImg) {
            bladeImg.src = bladeImage;
        } else {
            let _bladeImg = this.createElement({type: 'img', className: 'bladeImg'});
            _bladeImg.src = bladeImage;
            _bladeImg.style.pointerEvents = 'none'; // 设置刀片图不触发任何鼠标事件
            this.insertElement(this.editorBox, _bladeImg);
        }
        // 同步展示在预览区域
        let sampleData = {...sampleInfo};
        sampleData.src = sampleInfo.baseImgInfoList[0].baseImg;

        this.displayPreview('sample', sampleData, this.previewBox);
        this.displayPreview('sample', sampleData, this.exportObj);
    }

    // 根据class查找元素 @target 为要查找的元素父级
    findElementByClass(target, className) {
        if (!target) return;
        return this.searchNode(target, className);
    }

    hasEditImage() {
        return this.image;
    }

    //  删除元素
    delElement() {
        if (!this.onMoveObj) return;
        let objType = this.onMoveObj.getAttribute('data-type');
        let colorId = this.onMoveObj.getAttribute('colorId');
        // 删除色块
        if (objType === 'color') {
            this.delColorPiece(colorId);
            // 这里不能用强等，因为data里的是数字类型，元素的自定义属性取出来是字符串
            this.colorConfig.data = this.colorConfig.data.filter(item => item.colorId != colorId);
        }
        // 删除图片
        if (objType === 'pattern') {
            this.findElementByClass(this.elementsParent, 'editImage').remove();
            this.findElementByClass(this.editorBox, 'pattern-ctl').remove();
            this.findItemByAttr(this.findElementByClass(this.previewBox, 'preview-container'), 'pattern').remove();
        }
        if(objType === 'text'){
            this.delText();
        }
        // 删除后为了避免操作框还停留在页面上，所以要隐藏掉
        this.hideFrame();
        // 删除元素后置空当前操作的对象
        this.onMoveObj = null;
    }
    delText(){
        this.findElementByClass(this.elementsParent, 'edit-text').remove();
    }

    // 老版坐标点数转换
    posOldToNew(data) {
        let res = [];
        if (data) {
            data.forEach((item) => {
                let rs = item.split(',');
                res.push({pointx: Number(rs[0]) + 300, pointy: Number(rs[1]) + 300});
            })
        }
        return res;
    }

    // 根据转化后的老坐标获取贴图位置的大小和坐标
    getRectCenterOld(data) {
        const getPosition = (axis, type) => {
            return Math[type].apply(
                null,
                data.map((item) => item[axis])
            );
        };
        // top 最大值
        const mxTop = getPosition('pointy', 'max');
        // top最小值
        const minTop = getPosition('pointy', 'min');

        // left最大值
        const mxLeft = getPosition('pointx', 'max');
        // left最小值
        const minLeft = getPosition('pointx', 'min');

        // 根据坐标获取矩形宽高
        const height = mxTop - minTop - 10;
        const width = mxLeft - minLeft;
        return {left: minLeft, top: minTop, height, width};
    }

    // 合成图片,每次图片的编辑都会同步坐标合成单个图片
    async compositeImage({left, top, width, height, zIndex, colorId, transform}) {
        // 先检查页面有没有供合成的div；
        let canvas = document.getElementById('compositeCan');
        let imgForRender = this.findElementByClass(canvas, 'single-export');
        const setStyle = () => {
            let imgInRender = document.getElementsByClassName('single-export')[0];
            imgInRender.style.width = width + 'px';
            imgInRender.style.height = height + 'px';
            imgInRender.style.transform = transform;
            imgInRender.src = this.onMoveObj.getAttribute('src');
            imgInRender.style.zIndex = zIndex;
            imgInRender.style.left = left + 'px';
            imgInRender.style.position = 'absolute';
            imgInRender.style.top = top + 'px';
            imgInRender.style.pointerEvents = 'none';
        }
        if (!this.findElementByClass(canvas, 'single-export-container')) {
            let containerForRender = this.createElement({type: 'div', className: 'single-export-container'})
            containerForRender.style.position = 'relative';
            containerForRender.style.width = '100%';
            containerForRender.style.height = '100%';
            imgForRender = this.createElement({type: 'img', className: 'single-export'})
            this.insertElement(containerForRender, imgForRender);
            this.insertElement(canvas, containerForRender);
        }
        setStyle();

    }

    // 根据两个版本的坐标点返回坐标，和贴图区域大小
    getSliceSize(data,key='pathPoints',colorIndex) {
        let width, height, left, top;
        // 根据转后的坐标得到大小和位置
        if (data.newVersion) {
            let res = this.getRectCenter(this.blade.getCurrentBlade(this.currentBlade,colorIndex)[key]);
            width = res.width;
            height = res.height;
            left = res.left;
            top = res.top;
        } else {
            let old = this.posOldToNew(data.data);
            let res = this.getRectCenter(JSON.stringify(old));
            width = res.width;
            height = res.height;
            left = res.left;
            top = res.top;
        }
        return {
            width, height, left, top
        }
    }
    // 根据当前元素类型ID找预览区的元素
    findItemByTypeId(id){
        let type = this.getElementType();
        const elemnts = Array.from(this.previewImgContainer.childNodes);
        return elemnts.find(item => item.getAttribute(type+'id') == id);
    }

    // 根据当前元素类型ID找导出图的元素
    findItemByTypeIdInExport(id){
        let type = this.getElementType();
        const elemnts = Array.from(this.findElementByClass(this.exportObj, 'preview-container').childNodes);
        return elemnts.find(item => item.getAttribute(type+'id') == id);
    }
    // 合成图片,colorLen means the amount of img for one sample
    async doComposite(imgInfo, colorIndex,colorLen) {
        console.log(colorIndex,'colorIndexcolorIndexcolorIndexcolorIndex')
        const currentObjType = this.getElementType();
        if (currentObjType === 'pattern' || currentObjType === 'color') {
            // 导出图的DIV
            let exportParent = this.findElementByClass(document.getElementById('element-for-export'), 'preview-container');
            // 预览区的图片，因为图片只有一张，所以可以不用根据ID找
            let previewImg = this.findItemByAttr(this.previewImgContainer, currentObjType);
            // 根据当前操作对象的类型来找导出的DIV里面的所元素
            let exportImg = this.findItemByAttr(exportParent, currentObjType);
            // 如果当前预览区的图片没有ID 则根据 贴图区域的数量下标给值
            let layerId = currentObjType + colorIndex;

            // 如果是色块的预览区 要根据colorId找到当前操作的元素
            if(currentObjType === 'color'){
                previewImg = this.findItemByTypeId(this.onMoveObj.getAttribute('colorId'));
            }
            if(!previewImg.getAttribute('layerid')){
                previewImg.setAttribute('layerid',layerId);
            }
            // 如果是多刀片的话，可以插入多个图
            let isExist = Array.from(this.previewImgContainer.childNodes).find(item => item.getAttribute('layerid') == layerId);
            if (colorIndex > 0 && !isExist){
                previewImg = this.createElement({type: 'img'});
                previewImg.setAttribute('imgtype',currentObjType);
                previewImg.setAttribute('layerid', layerId);
                this.insertElement(this.previewImgContainer, previewImg);
            }
            if(isExist){
                previewImg = isExist;
            }

            // 这里是给导出效果图的DIV用的
            let isExistExport = Array.from(exportParent.childNodes).find(item => item.getAttribute('layerid') == layerId);
            if (colorIndex > 0 && !isExistExport){
                exportImg = this.createElement({type: 'img'});
                exportImg.setAttribute('imgtype',currentObjType);
                exportImg.setAttribute('layerid', layerId);
                this.insertElement(exportParent, exportImg);
            }
            if(currentObjType === 'color'){
                exportImg = this.findItemByTypeIdInExport(this.onMoveObj.getAttribute('colorId'));
            }
            // 导出图的DIV 也要这样操作一遍
            if(!exportImg.getAttribute('layerid')){
                exportImg.setAttribute('layerid',layerId);
            }

            if(isExistExport){
                exportImg = isExistExport;
            }
            let originImg = this.onMoveObj.getAttribute('src');

            // 获取当前操作对象的坐标和尺寸
            const {clientWidth, clientHeight, offsetLeft, offsetTop} = this.onMoveObj;

            // 数据格式转换，转换预览版和体验版的坐标
            let transferPos = this.getDesignSizeByPoints(colorIndex, 'pathPoints');// 先拿到坐标点的集合

            // 拿到切图的尺寸和坐标，切图和 贴图是两个不同的坐标
            const clipInfo = this.getSliceSize(transferPos,'clipPath',colorIndex);

            // 实例化透视以及角度算法的canvas
            const CanvasManager = new canvasManager({width: this.exportSize, height: this.exportSize})

            // 贴图坐标转换
            const canvasArea = stringToArray(transferPos.data, {width: this.exportSize, height: this.exportSize});

            // 创建一个canvas用来合成第一张图
            const {can, context} = this.createCan();
            let src, originSrc;
            // 设置画布的裁剪画布的大小来裁剪图片
            can.width = clipInfo.width;
            can.height = clipInfo.height;

            const imgForOrigin = new Image();
            imgForOrigin.src = originImg;
            context.translate(offsetLeft - (this.exportSize - can.width) /2 + (this.exportSize - this.canvasSize) / 2 ,  offsetTop - (this.exportSize - can.height) / 2 + + (this.exportSize - this.canvasSize) / 2)
            // 设置图片允许跨域
            imgForOrigin.setAttribute("crossOrigin", 'Anonymous');
            // 图片在画布的大小和位置

            imgForOrigin.onload = async () => {
                context.drawImage(imgForOrigin, 0, 0, clientWidth, clientHeight);
                originSrc = can.toDataURL('image/webp', 1);
                // 新数据的样板比老版本大一点点
                canvasArea.forEach((item) => {
                        item.x = transferPos.newVersion ? item.x * 1.3 : item.x * 0.9;
                        item.y = transferPos.newVersion ? item.y * 1.3 : item.y * 0.9;
                })
                previewImg.style.top = '0';
                previewImg.style.left = '0';
                previewImg.style.zIndex = imgInfo.zIndex;
                // 合成后的图赋值给预览图，和预备效果图合成的DIV
                previewImg.src = await CanvasManager.genImg(canvasArea, originSrc);
                previewImg.style.display = 'block';

                exportImg.src = previewImg.src;
                exportImg.style.display = 'block';
                exportImg.style.left = 0 + 'px';
                exportImg.style.top = 0 + 'px';
                exportImg.style.zIndex = imgInfo.zIndex;
                setTimeout(() => {
                    this.loadingPreview.style.display = 'none';
                },200)
            }
        }
    }

    // 创建canvas;
    createCan() {
        let can = document.createElement('canvas');
        let context = can.getContext('2d');
        return {
            can,
            context
        }
    }

    // 根据目标画布计算新比例
    getFitSizeByTarget(previewWidth) {
        let {offsetLeft, offsetTop, clientWidth, clientHeight} = this.onMoveObj; // 获取当前操作对象的x,y
        // 设计器 到预览区的坐标值转换，
        // 设计器的大小基准600， 预览区280 默认值
        const editorW = this.editorBox.offsetWidth;
        let moved = {
            left: '',
            top: '',
            width: '',
            height: '',
            zoom: '',
            zIndex: '',
            transform: '',
            innerText: '',
            textAlign: '',
            color: ''
        }
        // 计算出坐标值在编辑区的比例
        moved.left = offsetLeft / editorW * previewWidth + 'px';
        moved.top = offsetTop / editorW * previewWidth + 'px';
        moved.width = clientWidth / editorW * previewWidth + 'px';
        moved.height = clientHeight / editorW * previewWidth + 'px';
        moved.zIndex = this.onMoveObj.style.zIndex;
        moved.colorId = this.onMoveObj.getAttribute('colorid') || '';
        moved.transform = this.onMoveObj.style.transform;
        moved.innerText = this.onMoveObj.innerText;
        moved.textAlign = this.onMoveObj.style.textAlign;
        moved.color = this.onMoveObj.style.color;
        // moved.transform = 'rotateX(33deg) rotateY(17deg) rotateZ(349deg)';
        return moved;
    }
    // 找到预览区的指定元素 @type 为元素的类型，支持3种， color,patter,text
    findItemFromPreview(type){
        console.log()
        const currentObj = this.onMoveObj; // 当前操作的元素
        const id = currentObj.getAttribute(type+'id'); // 根据类型找到相应的ID
        return this.findItemByTypeId(id);
    }
    // 找到该类型元素的控制器
    findCtrFromEditBox(type, id){
        return Array.from(this.editorBox.childNodes).find(item => item.getAttribute(type+'id') === id)
    }
    // 元素操作后同步到预览区域 传入目标，比如预览区或者合成区
    asyncPosition(target) {
      // if(this.getElementType() ==='color') return;
        if (!this.onMoveObj) return;
        const previewBox = this.findElementByClass(target, 'preview-container');
        previewBox.style.perspective = '1200px';


        // 找到预览区的盒子，并且取盒子的大小
        const previewWidth = target.offsetWidth;
        let {offsetLeft, offsetTop, clientWidth, clientHeight} = this.onMoveObj; // 获取当前操作对象的x,y

        const targetType = this.getElementType(); // 判断当前操作的对象是图片还是色块

        // 根据当前元素类型找在预览区找到元素
        let _obj = this.findItemFromPreview(targetType);
        // console.log(_obj,'获取当前要转化坐标的对象', targetType, this.onMoveObj.getAttribute('colorid'))
        let {left, top, width, height, zIndex, colorId, transform} = this.getFitSizeByTarget(previewWidth);
        let currentObjId = this.onMoveObj.getAttribute(targetType +'id');
        let currentColorCtl = this.findCtrFromEditBox(targetType, currentObjId);
        console.log('找到控制器', currentColorCtl)
        let imgInfo = {
            zIndex: this.onMoveObj.style.zIndex,
            transform: this.onMoveObj.style.transform,
        }
        // 这里使用转换后的坐标和尺寸
        if (targetType === 'color' || targetType === 'pattern') {
            let ctl = targetType === 'color' ? 'color-ctl' : 'pattern-ctl';
            // 找到关联的控制框
            this.compositeImage(
                {
                    left: offsetLeft,
                    top: offsetTop,
                    width: clientWidth,
                    height: clientHeight,
                    zIndex: imgInfo.zIndex,
                    transform: imgInfo.transform
                }
            )
            // 调用合成方法
            setTimeout(async() => {
                this.defaultSample.baseImgInfoList[this.currentBlade].bladeInfoList.forEach((item,index) => {
                    this.doComposite(imgInfo,index,this.defaultSample.baseImgInfoList[this.currentBlade].bladeInfoList.length);
                })
            }, 100)
        }

        // 同步坐标和尺寸的时候要同时同步色块的操作框。
        if (targetType === 'color') {
            // 找到要替换颜色的色块ID
            let targetId = this.onMoveObj.getAttribute('colorId');
            let type = this.onMoveObj.getAttribute('colortype')
            // 删除前要保存当前色块的信息
            _obj.style.zIndex = imgInfo.zIndex;
            _obj.style.transform = imgInfo.transform;
            this.saveColorInfo(this.onMoveObj, type, targetId);
        }
        if(targetType === 'text'){
                this.asyncText(_obj);
        }
        if (!currentColorCtl) return;
        // 这里的坐标都是取当前操作对象的坐标
        currentColorCtl.style.width = clientWidth + 'px';
        currentColorCtl.style.height = clientHeight + 'px';
        currentColorCtl.style.left = this.onMoveObj.style.left;
        currentColorCtl.style.top = this.onMoveObj.style.top;
        currentColorCtl.style.transform = transform;
        this.asyncPositionForFrame(this.onMoveObj, offsetLeft, offsetTop);
        this.selectedInfo.info = this.createObjInfo();
        this.historyManage.add(this.onMoveObj);
    }

    asyncText(imgForPreview){
        const newPos = this.getFitSizeByTarget(280);
        imgForPreview.innerText = newPos.innerText;
        imgForPreview.style.textAlign = newPos.textAlign;
        imgForPreview.style.color = newPos.color;
        imgForPreview.style.left = newPos.left;
        imgForPreview.style.top = newPos.top;
        imgForPreview.style.zIndex = newPos.zIndex;
        imgForPreview.style.fontSize = newPos.fontSize;
        imgForPreview.style.position = 'absolute';
        this.loadingPreview.style.display = 'none';
    }

    // 如果没有选中的对象的时候显示画布的尺寸
    getCanvasSize() {
        this.selectedInfo.info = {
            width: 600,
            height: 600,
            zIndex: 0,
            left: 0,
            top: 0,
            type: '画布',
            exportSize: this.exportSize
        }
    }

    // 获取当前对象是色块还是图片
    getElementType() {
        return this.onMoveObj.getAttribute('data-type');
    }

    // 通过自定义属性找元素、
    findItemByAttr(parent, imageType) {
        return Array.from(parent.childNodes).find(item => item.getAttribute('imgtype') === imageType);
    }

    // 删除色块
    delColorPiece(id) {
        const deleteItem = (obj, target) => {
            Array.from(obj).forEach((item) => {
                // target === 'color-preview' && console.log(item)
                if (item.className === target && item.getAttribute('colorId') == id) {
                    item.remove();
                }
            });
        }

        let colorData = this.colorConfig.data.filter(item => item.colorId != id); // 删除当前色块
        let arr = [];
         // 去重
        colorData.forEach((item) => {
            if(arr.find(t => t.colorId != item.colorId)){
                arr.push(item);
            }
        })
        this.colorConfig.data = arr;
        deleteItem(this.editorBox.childNodes, 'color-ctl');// 色块的控制器
        deleteItem(this.elementsParent.childNodes, 'color-piece'); // 色块
        deleteItem(this.findElementByClass(this.previewBox, 'preview-container').childNodes, 'color-preview'); // 预览区的色块
    }

    // 生成色块的图片
    createColorPiece(colorType, color){
        let toDataURL = '';
        let canvas = document.createElement('canvas');
        let context = canvas.getContext('2d');
        canvas.height = this.colorConfig.size;
        canvas.width = this.colorConfig.size;
        context.fillStyle = color;
        if(colorType == 1) {
            context.arc(this.colorConfig.size/2, this.colorConfig.size/2, this.colorConfig.size/2, 0, Math.PI*2, true);
            context.fill();
        } else {
            context.fillRect(0, 0,this.colorConfig.size, this.colorConfig.size);
        }
        toDataURL = canvas.toDataURL('image/webp');
        return toDataURL;
    }
    // 设置色块的颜色
    setColor(color){
        // 找到要替换颜色的色块ID
        let targetId = this.onMoveObj.getAttribute('colorId');
        let type = this.onMoveObj.getAttribute('colortype');
        // 删除前要保存当前色块的信息
        this.saveColorInfo(this.onMoveObj, type, targetId);
        // 根据获取的ID获取到这个色块最后一次的尺寸和位置
        let lastColor = this.historyManage.find(this.onMoveObj);
        if(!lastColor){
            throw new Error('con not find colorId in function setColor. id:' + lastColor);
        }
        this.delColorPiece(targetId);
        this.addColorPiece(lastColor.type,color,lastColor);
    }
    // 查看有没有已经存在ID
    colorIdExist(_colorId){
        const nodes = Array.from(this.elementsParent.childNodes);
        nodes.forEach((item) => {
            if(item.getAttribute('colorId') == _colorId && item.className === 'color-piece'){
                item.remove();
            }
        })
    }
    // 显示预览区域的loading
    showLoading(){
        this.loadingPreview.style.display = 'block';
    }
    // 添加色块 type为色块类型，2为方形，1为圆形
    addColorPiece(type = 2, color = 'red',lastColor = null) {
        this.showLoading();
        // 如果是更换颜色需要沿用本身的type
        if(lastColor){
            type = lastColor.colorType;
        }
        if (!this.defaultSample) {
            alert('请选择样版');
            return;
        }
        if (this.colorConfig.data?.length === 5) {
            alert('最多只能添加五个色块！');
            return;
        }
        let colorImg = this.createColorPiece(type,color);
        const _editor = this.elementsParent;
        // 设置id,并且添加绑定事件，用于图层叠加可以选中色块
        let _colorId;
        if(this.colorConfig.data){
            _colorId = this.colorConfig.data.length + 1;
        }else{
            _colorId = 1;
        }
        // 如果是更换颜色需要沿用本身的id
        if(lastColor){
            _colorId = lastColor.id
        }
        // 如果有已经存在的ID，先删除掉
        this.colorIdExist();
        // 老规矩，没有找到就创建一个元素然后插入。
        let colorPiece = this.createElement({type: 'img', className: 'color-piece'});
        colorPiece.src = colorImg;

        // 替换的颜色的时候要把这个色块最后保存的信息赋值给新的色块
        if(lastColor){
            colorPiece.style.width =  lastColor.width + 'px';
            colorPiece.style.height = lastColor.width + 'px';
            colorPiece.style.left = lastColor.left + 'px';
            colorPiece.style.top = lastColor.top + 'px';
        }else{
            colorPiece.style.width =  this.colorConfig.size + 'px';
            colorPiece.style.height = this.colorConfig.size + 'px';
        }
        colorPiece.setAttribute('data-type', 'color');
        colorPiece.setAttribute('colorId', _colorId);

        colorPiece.setAttribute('colortype', type);
        colorPiece.draggable = false;
        this.insertElement(_editor,colorPiece)
        // 添加色块后选中的对象设置为色块并且设置操作框的大小和色块一致
        this.onMoveObj = colorPiece;
        this.frame.style.width = this.colorConfig.size + 'px';
        this.frame.style.height = this.colorConfig.size + 'px';
        colorPiece.style.zIndex = this.layers?.length + 1;
        this.saveColorInfo(colorPiece, type, _colorId);

        // 添加完色块，在操作框同级插入隐藏的边框用于图层叠加选不中下一层用。
        let colorCtl = this.createElement({
            type: 'div',
            className: 'color-ctl',
        });

        if(lastColor){
            colorPiece.style.width =  lastColor.width + 'px';
            colorPiece.style.height = lastColor.width + 'px';
            colorPiece.style.left = lastColor.left + 'px';
            colorPiece.style.top = lastColor.top + 'px';
        }else{
            colorCtl.style.width = this.colorConfig.size + 'px';
            colorCtl.style.height = this.colorConfig.size + 'px';
        }

        colorCtl.style.zIndex = '12';
        colorCtl.style.pointerEvents = 'auto';
        colorCtl.setAttribute('el-type', 'ctl');

        colorCtl.setAttribute('colorId', _colorId);
        this.layers.push(colorPiece.style.zIndex); // 保存起来作为管理层级用
        // 这个方法用于色块和色块叠加的时候 选不中色块的特殊处理。
        colorCtl.addEventListener('mousedown', () => {
            setTimeout(() => {
                colorCtl.style.display = 'none';
                this.onMoveObj = colorPiece;
            })
        })
        this.insertElement(this.editorBox, colorCtl)
        // this.addEventForElement(colorPiece)
        this.resetFramePosition(0, 0); // 重置操作框的坐标
        this.displayPreview('color', colorPiece, this.previewBox);
        this.displayPreview('color', colorPiece, this.exportObj);
    }

    // 选中图片的时候需要显示所有色块关联的操作框。
    showAllColorCtl(state) {
        setTimeout(() => {
            Array.from(this.editorBox.childNodes).forEach((item, index) => {
                if (item.className === 'color-ctl') {
                    // console.log(index);
                    item.style.display = state;
                    item.style.left = this.colorConfig.data[item.getAttribute('colorId') - 1]?.left;
                    item.style.top = this.colorConfig.data[item.getAttribute('colorId') - 1]?.top;
                }
            });

            // 图片控制器
            const imageCtl = this.findElementByClass(this.editorBox, 'pattern-ctl');
            imageCtl && (imageCtl.style.display = state);
        }, 500)
    }

    // 保存色块信息
    saveColorInfo(colorPiece, type,_colorId) {
        this.colorConfig.data.push({
            colorId: _colorId,
            color: colorPiece.style.backgroundColor,
            size: {
                h: colorPiece.clientHeight,
                w: colorPiece.clientWidth,
            },
            type: type,
            left: colorPiece.offsetLeft +'px',
            top: colorPiece.offsetTop + 'px'
        })
    }
    // 设置颜色
    setTextColor(val){
        this.onMoveObj.style.color = val;
        const imageInfo = {};
        this.displayPreview('text', imageInfo, this.previewBox);
        this.displayPreview('text', imageInfo, this.exportObj);
    }
    // 设置文字对齐方式
    setAlign(val){
        this.onMoveObj.style.textAlign = val;
        const imageInfo = {};
        this.displayPreview('text', imageInfo, this.previewBox);
        this.displayPreview('text', imageInfo, this.exportObj);
    }
    // 添加文字
    addText(){
        let elementForText = this.createElement({type:'div', className:'edit-text'});
        elementForText.innerText = '请输入文字';
        elementForText.contentEditable = 'true';
        elementForText.style.left = 600 / 2 + 'px';
        elementForText.style.top = 600 / 2 + 'px';
        elementForText.style.zIndex = 600 / 2;
        elementForText.setAttribute('textid', 1);
        elementForText.setAttribute('data-type','text');
        const imageInfo = {};
        this.insertElement(this.elementsParent,elementForText);

        // 添加完文字，在操作框同级插入隐藏的边框用于图层叠加选不中下一层用。
        let textCtl = this.createElement({
            type: 'div',
            className: 'text-ctl',
        });
        textCtl.style.zIndex = '20';
        textCtl.style.pointerEvents = 'auto';
        textCtl.setAttribute('textid', 1);
        this.onMoveObj = elementForText;
        this.insertElement(this.editorBox, textCtl)
        this.displayPreview('text', imageInfo, this.previewBox);
        this.displayPreview('text', imageInfo, this.exportObj);
    }
    // 添加图片
    addImage(imageInfo) {
        // 这里要处理一下，不然找不到src
        if (!imageInfo) {
            alert('请传入有效的图片信息');
            return;
        }
        if (!this.defaultSample) {
            alert('请选择样版');
            return;
        }
        this.showLoading();
        this.designContents.bladeList[this.currentBlade].contents.images = {...imageInfo};
        // 如果是第一次添加图案，需要生成img标签再插入到父节点，
        let currentImage = this.findElementByClass(this.elementsParent, 'editImage');

        // 先判断设计内容里面是否已经有图片了，如果有的话就替换掉
        if (!currentImage) {
            currentImage = this.createElement({type: 'img', className: 'editImage'})
            currentImage.src = imageInfo.canvasImg || imageInfo.src;
            this.insertElement(this.elementsParent, currentImage);
            this.image = currentImage;
            this.onMoveObj = currentImage;
            currentImage.draggable = false;
            currentImage.setAttribute('patternid',1);
            currentImage.setAttribute('data-type','pattern');
            currentImage.style.zIndex = 1;
            this.layers.push(currentImage.style.zIndex); // 保存起来作为管理层级用
        } else {
            // 更换图片需要重置图片对象的宽度。不然会压缩变形
            currentImage.src = imageInfo.canvasImg || imageInfo.src;
            currentImage.style.width = '100%';
            currentImage.style.height = 'auto';
            this.findElementByClass(this.editorBox, 'pattern-ctl').remove();
        }
        // 添加图片，在操作框同级插入隐藏的边框用于图层叠加选不中下一层用。
        let imageCtl = this.createElement({
            type: 'div',
            className: 'pattern-ctl',
        });
        imageCtl.style.width = currentImage.style.width;
        imageCtl.style.height = currentImage.style.height;
        imageCtl.style.zIndex = '12';
        imageCtl.style.pointerEvents = 'auto';
        imageCtl.setAttribute('el-type', 'ctl');
        imageCtl.setAttribute('patternid', '1');

        // 这个方法用于色块和色块叠加的时候 选不中色块的特殊处理。
        imageCtl.addEventListener('mousedown', () => {
            setTimeout(() => {
                imageCtl.style.display = 'none';
                this.onMoveObj = currentImage;
            })
        })
        this.insertElement(this.editorBox, imageCtl);
        currentImage.onload = () => {
            this.setImageToFit(currentImage);
            // 因为要同时把设计器内的编辑内容同步到 预合成DIV和预览区，所以要调用两次
            this.displayPreview('pattern', imageInfo, this.previewBox);
            this.displayPreview('pattern', imageInfo, this.exportObj);
            this.asyncPosition(this.previewBox);
        }
    }
    // 聚焦到点击的元素
    focusSelectedItem(obj) {
        if (obj.id === 'move') {
            return;
        }
        // 设置当前操作的对象
        this.onMoveObj = obj;
        // 同步当前操作对象的尺寸以及坐标到操作框
        this.frame.style.width = obj.clientWidth - 4 + 'px';
        this.frame.style.height = obj.clientHeight - 2 + 'px';
        this.frame.style.left = obj.offsetLeft + 'px';
        this.frame.style.top = obj.offsetTop + 'px';
    }

    // 创建当前对象的信息
    createObjInfo() {
        console.log(this.getElementType(),'type')
        const {clientWidth, clientHeight, offsetLeft, offsetTop} = this.onMoveObj;
        const _this = this;
        const info = {
            type: _this.getElementType(),
            width: clientWidth,
            height: clientHeight,
            left: offsetLeft,
            top: offsetTop,
            zIndex: _this.onMoveObj.style.zIndex,
            exportSize: this.exportSize
        }
        return info;
    }

    // 设置图片居中
    setImageToCenter(target, width, height) {
        let _left = target / 2 - width / 2,
            _top = target / 2 - height / 2;
        this.onMoveObj.style.left = _left + 'px';
        this.onMoveObj.style.top = _top + 'px';
    }

    // 设置不同尺寸的图片适应到画布的尺寸，比如超长图，要按照比例缩放
    setImageToFit(currentImage) {
        const {clientWidth, clientHeight} = currentImage;
        const basicSize = 600;
        // 因为画布大小是固定600
        const imageOverSize = [clientWidth, clientHeight];
        let zoom = 1; // 默认缩放比例为100%;
        // 判断图片的宽高是否有一个值大于画布
        const overSize = imageOverSize.find(item => item > basicSize);
        let _left, _top;

        // 图片原始尺寸超出画布
        if (overSize) {
            zoom = zoom - (overSize - basicSize) / basicSize; // 根据超出的数值计算应该缩小多少比例
            // 得出适应画布的宽高
            let fitWidth = clientWidth * zoom;
            let fitHeight = clientHeight * zoom;
            // 设置宽高和坐标
            currentImage.style.width = fitWidth + 'px';
            currentImage.style.height = fitHeight + 'px';
            this.setImageToCenter(basicSize, fitWidth, fitHeight)
        } else {
            // 图片原始大小小于或者等于画布大小
            _left = basicSize / 2 - clientWidth / 2;
            _top = basicSize / 2 - clientHeight / 2;
            // 设置图片居中
            currentImage.style.left = _left + 'px';
            currentImage.style.top = _top + 'px';
        }

        this.asyncPositionForFrame(currentImage, _left, _top);
    }

    // 获取当前需要的坐标点数
    getDesignSizeByPoints(index = 0,key='pathPoints') {
        let res = {newVersion: false, data: []};
        if (this.defaultSample) {
            let currentPoints = this.blade.getCurrentBlade(this.currentBlade,index);
            if (currentPoints.pathPoints) { // 判断是不是新坐标标点的
                res.newVersion = true;
                res.data = JSON.parse(currentPoints[key]);
                res.data = res.data.map(element => {
                    element.pointx = parseInt(element.pointx);
                    element.pointy = parseInt(element.pointy);
                    element.pointx = element.pointx - this.canvasSize /2;
                    element.pointy = element.pointy - this.canvasSize /2;
                    return element.pointx + ',' + element.pointy;
                })
                res.data.pop();
            } else {
                // 旧数据
                res.data = [currentPoints.p1, currentPoints.p2, currentPoints.p3, currentPoints.p4];
            }
            return res;
        }
        return [];
    }

    // 计算每张刀片裁剪位置
    getRectCenter(path) {
        if (!path) {
            return {left: 0, top: 0, height: 600, width: 600}
        }
        // 先格式化一下，因为拿到的时候是字符串
        const formatPath = (path) => {
            path = eval(path);
            let formatedPath = [];
            path.forEach((item) => {
                formatedPath.push({
                    x: item.pointx,
                    y: item.pointy,
                });
            });
            return formatedPath;
        }
        if (typeof path === 'string') {
            path = formatPath(path);
        }
        /*
         * axis 为坐标轴
         * type min 为最小值
         * type max 为最大值
         * */
        const getPosition = (axis, type) => {
            return Math[type].apply(
                null,
                path.map((item) => item[axis])
            );
        };
        // top 最大值
        const mxTop = getPosition('y', 'max');
        // top最小值
        const minTop = getPosition('y', 'min');

        // left最大值
        const mxLeft = getPosition('x', 'max');
        // left最小值
        const minLeft = getPosition('x', 'min');

        // 根据坐标获取矩形宽高
        const height = mxTop - minTop;
        const width = mxLeft - minLeft;

        // 计算矩形中心
        const centerY = height / 2 + minTop;
        const centerX = width / 2 + minLeft;
        return {left: minLeft, top: minTop, height, width};
    }

    // 添加图片的时候需要让外框大小和图片同步
    asyncPositionForFrame(currentImage, left, top) {
        const {clientWidth, clientHeight} = currentImage;
        this.frame.style.width = clientWidth - 4 + 'px'; // 这里减去的数值 是因为有一点点偏差
        this.frame.style.height = clientHeight - 2 + 'px';
        this.frame.style.transform = currentImage.style.transform;
        this.showFrame();
        // 添加或者替换图片后，图片和操作框的位置都要恢复到默认
        this.resetFramePosition(left, top);
    }

    // 渲染预览区域
    displayPreview(type, element, areaTarget) {
        let previewImgContainer = this.searchNode(areaTarget, 'preview-container');
        // 先查找是否包含了preview-container
        if (!previewImgContainer) {
            // 如果没有则创建一个
            previewImgContainer = this.createElement({type: 'div', className: 'preview-container'});
            this.insertElement(areaTarget, previewImgContainer);
            this.previewImgContainer = previewImgContainer;
        }

        // 找到预览区域的元素
        const findItemInPreview = (el) => {
            let _container = this.searchNode(areaTarget, 'preview-container');
            return Array.from(_container.childNodes).find(item => item.getAttribute('imgType') === el);
        }
        // 给预览区创建元素
        const createElementForPreview = (type) => {
            let tagType = type === 'text' ? 'div' : 'img';
            let imgForPreview = this.createElement({type:  tagType});
            // 设置自定义属性来区分图案和色块
            imgForPreview.setAttribute('imgType', type);
            if (type === 'color') {
                imgForPreview.className = 'color-preview';
            }
            if (type === 'sample') {
                imgForPreview.style.zIndex = '20'; // 因为样板总是要在最上面一层
            }
            if(type !== 'sample'){
                imgForPreview.setAttribute(type+'id',this.onMoveObj.getAttribute(type+'id'));
            }
            if(type === 'text'){
                 this.asyncText(imgForPreview);
            }
            imgForPreview.draggable = false
            imgForPreview.src = element.canvasImg || element.src;
            this.insertElement(previewImgContainer, imgForPreview); // 在预览区域插入图片
        }
        switch (type) {
            case 'pattern':
                // 插入图片需要先判断如果有图片直接替换图片路径即可，因为业务目前一个刀片只允许存在一张图片
                let _pattern = findItemInPreview(type);
                if (_pattern) {
                    _pattern.src = element.canvasImg || element.src;
                } else {
                    createElementForPreview(type);
                }
                break;
            case 'sample':
                let _sample = findItemInPreview(type);
                if (_sample) {
                    _sample.src = element.canvasImg || element.src;
                } else {
                    createElementForPreview(type);
                }

                // 用于效果图的DIV
                const exportElements = Array.from(this.findElementByClass(this.exportObj, 'preview-container').childNodes);
                let exportSample = exportElements.find(item => item.getAttribute('imgType') == type);
                if (exportSample) {
                    exportSample.src = element.canvasImg || element.src;
                } else {
                    this.createSampleForExport(element.canvasImg || element.src);
                }
                break;
            case 'color':
                createElementForPreview(type);
                break;
            case 'text':
                createElementForPreview(type);
                break;
        }
        // 更换图片后需要同步下尺寸和大小到 预览区
        this.asyncPosition(areaTarget);
    }

    // 给导出效果图的DIV添加样板
    createSampleForExport(src) {
        let _element = this.createElement({type: 'img',});
        _element.draggable = false;
        _element.setAttribute('imgtype', 'sample');
        _element.src = src
        this.insertElement(this.exportObj, _element);
    }

    // 根据className查找节点是否存在
    searchNode(target, node) {
        return Array.from(target.childNodes).find(item => item.className === node);
    }

    // 创建外框控制的父元素 并且添加8个操作点
    createFramePoints() {
        const _this = this;
        let editorCtrl = document.createElement("div");
        editorCtrl.className = this.frameParentName;
        const parentEl = document.getElementById('move');
        // 生成编辑框的操作点并插入到控制对象的元素中
        this.framePoints.forEach((item) => {
            let point = document.createElement("span");
            point.className = item;
            const setDefaultPosition = (e) => {
                _this.pointA.x = e.pageX;
                _this.pointA.y = e.pageY;
            }
            switch (item) {
                case 'top-right':
                    point.addEventListener('click', () => {
                        _this.delElement();
                    })
                    break;
                case 'top':
                    point.onmousedown = function (e) {
                        _this.isTopHold = true;
                        setDefaultPosition(e);
                    }
                    break;
                case 'left':
                    point.onmousedown = function (e) {
                        _this.isLeftHold = true;
                        console.log(e.pageX, e.pageY);
                        setDefaultPosition(e);
                    }
                    break;
                case 'right':
                    point.onmousedown = function (e) {
                        _this.isRightHold = true;
                        setDefaultPosition(e);
                    }
                    break;
                case 'bottom':
                    point.onmousedown = function (e) {
                        _this.isBottomHold = true;
                        setDefaultPosition(e);
                    }
                    break;
            }

            editorCtrl.appendChild(point);
        })
        parentEl.appendChild(editorCtrl);
        this.createRotationBtn(parentEl);
        this.addZoomListener();
    }

    // 创建旋转按钮
    createRotationBtn(container) {
        const _this = this;
        const rotationCtl = this.createElement({type: 'span', id: 'rotationCtl'});

        rotationCtl.onmousedown = function (e) {
            _this.currentAction = 'rotate';
            _this.isHold = true;
            if (_this.allAngel === 0) {
                _this.pointA.x = _this.editorBox.clientWidth / 2 + _this.editorBox.offsetLeft;
                _this.pointA.y = _this.editorBox.clientHeight / 2 + _this.editorBox.offsetTop;
            }
            // 获取鼠标的起始点坐标
            if (_this.count < 1) {
                _this.pointB.x = e.pageX;
                _this.pointB.y = e.pageY;
                _this.count++;
            }
        }
        container.appendChild(rotationCtl);
    }

    // 设置操作框到默认位置
    resetFramePosition(left, top) {
        // 这里的100 和 20 是固定要减去的值。和鼠标事件里的算法相同
        this.frame.style.left = left + 'px';
        this.frame.style.top = top + 'px';
    }

    // 添加图片放大缩小监听
    addZoomListener() {
        this.editorBox.addEventListener("mousewheel", (event) => {
            let delta = 0;
            if (!event) event = window.event; // 为了兼容浏览器
            if (event.wheelDelta) {
                delta = event.wheelDelta / 120;
                if (window.opera) delta = -delta;
            } else if (event.detail) {
                delta = -event.detail / 3;
            }
            let width = this.onMoveObj.offsetWidth;
            let height = this.onMoveObj.offsetHeight;
            let zoomInWidth = width * 1.1 + "px";
            let zoomInHeight = height * 1.1 + "px";
            let zoomOutWidth = width * 0.9 + "px";
            let zoomOutHeight = height * 0.9 + "px";
            if (delta > 0) {
                this.onMoveObj.style.width = zoomInWidth;
                this.onMoveObj.style.height = zoomInHeight;
                this.frame.style.width = zoomInWidth;
                this.frame.style.height = zoomInHeight;
            } else if (delta < 0) {
                if (width > 200) {
                    this.onMoveObj.style.width = zoomOutWidth;
                    this.onMoveObj.style.height = zoomOutHeight;
                    this.frame.style.width = zoomOutWidth;
                    this.frame.style.height = zoomOutHeight;
                }
            }
            this.asyncPosition(this.previewBox); // 同步到预览区
            // this.asyncPosition(this.exportObj); // 同步到合成区
            // this.asyncPositionForFrame() // 同步外边框
        }, false);

    }

    // 清除编辑框的移动监听事件
    removeFrameMouseEvent() {
        this.frame.onmousedown = null;
        this.frame.style.pointerEvents = 'none';
        this.frame.onmousemove = null;
        this.frame.onmouseup = null;
    }

    // 导出效果图
    async exportImage() {
        if(!domtoimage){
            throw new Error('domtoimage is not installed!');
        }
        if (!this.onMoveObj) return false;
        const nodeElement = this.findElementByClass(this.exportObj, 'preview-container');
        let base64Img = await domtoimage.toPng(nodeElement,{
            width: this.exportSize,
            height: this.exportSize,
        });
        return base64Img
    }
    // gen images for toolKit actions
    genImgForAction(action){
        // 找到合成单图的元素
        let singleComposite = document.getElementById('compositeCan');
        let singleCompBox = this.findElementByClass(singleComposite,'single-export-container');
        let singleCompImg = this.findElementByClass(singleCompBox,'single-export');

        // 找到预览区域的图获取SRC赋值给单图片合成
        let previewImg = this.findItemByAttr(this.previewImgContainer, 'pattern');
        singleCompImg.src = previewImg.src;
        const transform = this.onMoveObj.style.transform;
        console.log('gen images for toolKit actions',transform)
        if(!transform) return;
        // 找到最终导出效果的div
        let exportDiv = this.findElementByClass(this.exportObj,'preview-container');
        let exportImg = this.findItemByAttr(exportDiv, 'pattern');
        exportImg.style.transform = transform;

        // 把编辑区操后的旋转角度给合成单图的图片
        previewImg.style.transform = transform;
        singleCompImg.style.transform = transform;
    }
    // 工具集合
    toolKit(key) {
        if (!this.onMoveObj) return
        const canvasWidth = 600; // 编辑器画布大小为固定600
        const currentObj = this.onMoveObj;
        let {clientWidth, clientHeight} = currentObj;
        let currentIndex = parseInt(currentObj.style.zIndex);
        const currentObjType = this.getElementType(); // 获取当前操作是图片还是色块
        let currentColorCtl = null;
        // 根据当前对象的类型来调整层级
        if (currentObjType === 'image') {
            currentColorCtl = Array.from(this.editorBox.childNodes).find(item => item.className === 'pattern-ctl');
        }
        if (currentObjType === 'color') {
            currentColorCtl = Array.from(this.editorBox.childNodes).find(item => item.getAttribute('colorid') === currentObj.getAttribute('colorid'));
        }
        const actions = {
            // 上一层
            layerUp: () => {
                // 编辑内容部分
                const designElements = Array.from(this.elementsParent.childNodes);
                let currentEl = designElements.find(item => item.style.zIndex == currentIndex);
                console.log(currentEl, '当前的编辑内容', currentIndex)
                // 循环编辑内容的元素 找到比自己层级大，当前的元素的层级设置为比他层级大的那一层的值 + 1即可
                for (let i = 0; i < designElements.length; i++) {
                    if (designElements[i]?.style.zIndex >= currentEl.style.zIndex) {
                        currentEl.style.zIndex = parseInt(designElements[i].style.zIndex) + 1;
                        currentColorCtl && (currentColorCtl.style.zIndex = parseInt(currentEl.style.zIndex) + 10);
                        break;
                    }
                }
            },
            // 下一层
            layerDown: () => {
                const designElements = Array.from(this.elementsParent.childNodes);
                let currentEl = designElements.find(item => item.style.zIndex == currentIndex);
                console.log(currentEl, '当前的编辑内容', currentIndex)
                // 循环编辑内容的元素 找到比自己层级大，当前的元素的层级设置为比他层级大的那一层的值 + 1即可
                for (let i = 0; i < designElements.length; i++) {
                    if (designElements[i]?.style.zIndex <= currentEl.style.zIndex) {
                        currentEl.style.zIndex = parseInt(designElements[i].style.zIndex) - 1;
                        currentColorCtl && (currentColorCtl.style.zIndex = parseInt(currentEl.style.zIndex) + 10);
                        break;
                    }
                }
            },
            // 铺满区域
            fullCover: () => {
                // 先获取要操作对象的高和宽
                // 先找到比画布大的值
                let smallerVal = [clientWidth, clientHeight].find(item => item < canvasWidth);

                // 计算出差值的百分比
                let diff = canvasWidth / smallerVal;
                // 按照比例放大图片到铺满效果
                currentObj.style.width = clientWidth * diff + 'px';
                currentObj.style.height = clientHeight * diff + 'px';

                // 调用方法设置对象居中
                this.setImageToCenter(canvasWidth, clientWidth * diff, clientHeight * diff);
                // this.asyncPositionForFrame();

            },
            // 撑满宽度
            maxWidth: () => {
                currentObj.style.width = canvasWidth + 'px';
                this.setImageToCenter(canvasWidth, currentObj.clientWidth, currentObj.clientHeight);
            },
            // 撑满高度
            maxHeight: () => {
                currentObj.style.height = canvasWidth + 'px';
                this.setImageToCenter(canvasWidth, currentObj.clientWidth, currentObj.clientHeight);
            },
            // 水平镜像
            mirrorX: async () => {
                currentObj.style.transform = currentObj.style.transform === 'rotateY(180deg)' ? 'rotateY(0deg)' : 'rotateY(180deg)';
                currentObj.src = await domtoimage.toPng(currentObj,{
                            quality: '1',
                            width: this.canvasSize,
                            height: this.canvasSize,
                        });
                currentObj.style.transform = '';
            },
            mirrorY: async () => {
                currentObj.style.transform = currentObj.style.transform === 'rotateX(180deg)' ? 'rotateX(0deg)' : 'rotateX(180deg)';
                currentObj.src = await domtoimage.toPng(currentObj,{
                    quality: '1',
                    width: this.canvasSize,
                    height: this.canvasSize,
                });
                currentObj.style.transform = '';
            },
            // 水平居中
            centerW: () => {
                let _left = canvasWidth / 2 - currentObj.clientWidth / 2;
                currentObj.style.left = _left + 'px';
            },
            // 垂直居中
            centerH: () => {
                let _top = canvasWidth / 2 - currentObj.clientHeight / 2;
                currentObj.style.top = _top + 'px';
            },
            // 图案重复循环
            repeat: async () => {
                const currentObjType = this.getElementType();
                if (currentObjType === 'color') {
                    alert('图案不支持循环操作')
                    return;
                }
                // 创建一个隐藏的DIV 用于图案循环完毕后直接用canvas导出成图片替换掉原来的图片
                const boxForRepeat = this.createElement({type: 'div', className: 'repeat-box'});
                let img = this.onMoveObj.src;
                boxForRepeat.style.width = '600px'; // 这里是编辑区的大小，这里是写死的。
                boxForRepeat.style.height = '600px';

                // 把图片设置为背景，根据当前操作对象的宽度进行图案循环
                boxForRepeat.style.background = `url(${img}) repeat`;
                boxForRepeat.style.backgroundSize = this.onMoveObj.clientWidth + 'px';

                // 插入对象到页面并且导出成canvas图 替换掉设计区域的原图
                document.body.appendChild(boxForRepeat);
                let canvasImg = await domtoimage.toPng(boxForRepeat);
                Object.assign(this.designContents.bladeList[this.currentBlade].contents.images, {canvasImg: canvasImg});
                this.addImage(this.designContents.bladeList[this.currentBlade].contents.images);
                boxForRepeat.remove(); // 最后销毁掉
            }

        }
        actions[key] && actions[key]();
        if (key !== 'exportImage') {
            this.showLoading()
            this.asyncPosition(this.previewBox);
        }
        // this.asyncPosition(this.exportObj);
    }

    // 旋转事件
    onRotate(e) {
        const _this = this;
        _this.pointC.x = e.pageX;
        _this.pointC.y = e.pageY; // 获取结束点坐标
        // 计算旋转的角度
        let AB = {}, AC = {};

        AB.x = (_this.pointB.x - _this.pointA.x);
        AB.y = (_this.pointB.y - _this.pointA.y);
        AC.x = (_this.pointC.x - _this.pointA.x);
        AC.y = (_this.pointC.y - _this.pointA.y);
        // AC和AB 求出逆时针，还是顺时针
        let direct = (AB.x * AC.y) - (AB.y * AC.x);
        let lengthAB = Math.sqrt(Math.pow(_this.pointA.x - _this.pointB.x, 2) +
            Math.pow(_this.pointA.y - _this.pointB.y, 2)),
            lengthAC = Math.sqrt(Math.pow(_this.pointA.x - _this.pointC.x, 2) +
                Math.pow(_this.pointA.y - _this.pointC.y, 2)),
            lengthBC = Math.sqrt(Math.pow(_this.pointB.x - _this.pointC.x, 2) +
                Math.pow(_this.pointB.y - _this.pointC.y, 2));
        let cosA = (Math.pow(lengthAB, 2) + Math.pow(lengthAC, 2) - Math.pow(lengthBC, 2)) /
            (2 * lengthAB * lengthAC); // 余弦定理秋出旋转角
        let angleA = Math.round(Math.acos(cosA) * 180 / Math.PI);
        let allAngel = 0;
        if (direct < 0) {
            allAngel = -angleA; // 为负数的时候代表逆时针旋转
        } else {
            allAngel = angleA; // 正数为顺时针
        }

        // 如果之前旋转过要加上之前的旋转角度
        if (_this.allAngel) {
            allAngel += _this.allAngel;
        }
        _this.onMoveObj.style.transform = `rotate(${allAngel}deg)`;
        _this.frame.style.transform = `rotate(${allAngel}deg)`;
        _this.onMoveObj.setAttribute('rotate', allAngel);
        // _this.allAngel = allAngel;
        // _this.asyncPosition(_this.previewBox);
    }

    // 移动事件
    onMove(e) {
        const _this = this;
        e = e || window.event; // 兼容IE写法
        _this.showAllColorCtl('none');

        // 获取鼠标位置
        let x = e.clientX,
            y = e.clientY;
        // console.log('当前拖拽的目标',_this.onMoveObj)
        // 修改图片位置
        if (!_this.onMoveObj) return;
        _this.onMoveObj.style.left = x - _this.disX + "px";
        _this.onMoveObj.style.top = y - _this.disY + "px";

        // 记录最后的位置保存起来
        let currentColorId = _this.onMoveObj.getAttribute('colorId');
        if (currentColorId) {
            let currentIndex = currentColorId - 1;
            if (_this.colorConfig.data[currentIndex]) {
                _this.colorConfig.data[currentIndex].left = x - _this.disX + "px";
                _this.colorConfig.data[currentIndex].top = y - _this.disY + "px";
            }
        } else {
            let type = _this.getElementType();
            let currentImageCtl = _this.findElementByClass(_this.editorBox, type+'-ctl');
            currentImageCtl.style.left = x - _this.disX + "px";
            currentImageCtl.style.top = y - _this.disY + "px";
        }
        _this.frame.style.left = _this.defaultLeft + x - _this.disX + "px";
        _this.frame.style.top = _this.defaultTop + y - _this.disY + "px";
        // _this.frame.style.marginLeft = "0";
    }

    // 操作框点鼠标拖动事件
    topOnMove(e) {
        // 获取鼠标按下的时候鼠标的坐标
        let dis = {...this.pointA};

        // 获取当前鼠标停留的位置
        let lastY = e.pageY;

        // 新高度 = 原始高度 + 鼠标初始点 - 最后停留的位置
        if (this.objDefParams.height + (dis.y - lastY) < 10) return; // 边界处理
        this.onMoveObj.style.height = this.objDefParams.height + (dis.y - lastY) + 'px';

        // y坐标 = 原始y - 鼠标初始点 - 最后停留的位置
        this.onMoveObj.style.top = this.objDefParams.top - (dis.y - lastY) + 'px';
        this.asyncPositionForFrame(this.onMoveObj, this.onMoveObj.offsetLeft, this.onMoveObj.offsetTop);
        // this.asyncPosition(this.previewBox);

    }

    // 操作框点鼠标拖动事件
    leftOnMove(e) {
        // 获取当前鼠标停留的位置
        let lastX = e.pageX;

        // 第一次移动的位置
        let dis = {...this.pointA};
        let width, left;
        let diff = lastX - dis.x; // 计算第一次的A点到B点的差距
        // 鼠标往右移动代表缩小，差值会为正数
        if (diff > 0) {
            width = this.objDefParams.width - diff + 'px';
            left = this.objDefParams.left + diff + 'px';
        } else {
            // 鼠标往左数值会出现负数，则需要倒过来用oldX - lastX
            width = this.objDefParams.width + (dis.x - lastX) + 'px';
            left = this.objDefParams.left - (dis.x - lastX) + 'px';
        }
        if (width.split('px')[0] < 10) return; // 边界处理
        this.onMoveObj.style.width = width;
        this.onMoveObj.style.left = left;
        this.asyncPositionForFrame(this.onMoveObj, this.onMoveObj.offsetLeft, this.onMoveObj.offsetTop);
        // this.asyncPosition(this.previewBox);
    }

    // 操作框点鼠标拖动事件
    rightOnMove(e) {
        // 获取当前鼠标停留的位置
        let lastX = e.pageX;

        // 第一次移动的位置
        let dis = {...this.pointA};
        let width, left;
        let diff = lastX - dis.x; // 计算第一次的A点到B点的差距


        // 右边控制点 往右拉动不用计算left，因为图片本身默认是从左往右伸展的
        if (diff > 0) {
            width = this.objDefParams.width + diff + 'px';
        } else {
            width = this.objDefParams.width - (dis.x - lastX) + 'px';
        }
        if (width.split('px')[0] < 10) return; // 边界处理
        // 计算出新的宽度和坐标，因为要固定X轴的位置，所以x轴也要做计算
        this.onMoveObj.style.width = width;
        this.asyncPositionForFrame(this.onMoveObj, this.onMoveObj.offsetLeft, this.onMoveObj.offsetTop);
        // this.asyncPosition(this.previewBox);
    }

    // 操作框点鼠标拖动事件
    bottomOnMove(e) {
        // 获取鼠标按下的时候鼠标的坐标
        let dis = {...this.pointA};

        // 获取当前鼠标停留的位置
        let lastY = e.pageY;

        // 新高度 = 原始高度  + 鼠标初始点 - 最后停留的位置
        if (this.objDefParams.height + (lastY - dis.y) < 10) return; // 边界处理
        this.onMoveObj.style.height = this.objDefParams.height + (lastY - dis.y) + 'px';
        // y坐标 = 原始y - 鼠标初始点 - 最后停留的位置
        this.asyncPositionForFrame(this.onMoveObj, this.onMoveObj.offsetLeft, this.onMoveObj.offsetTop);
        // this.asyncPosition(this.previewBox);
    }

}