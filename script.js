let height = +document.getElementById('matrix-height').value;
let width = +document.getElementById('matrix-width').value;
let nodes = [];
let brs = [];
let matrix = [];
let ref = [];
let rref = [];
let kernel = [];
let pivotColumns = [];
let transposed = [];

let resultNodes = [];

function updateInput() {
    //////////////////////////////////////////////////////////
    for (let i = 0; i < nodes.length; i++) {
        for (let j = 0; j < nodes[i].length; j++) {
            document.body.removeChild(nodes[i][j]);
        }
    }
    for (let i = 0; i < brs.length; i++) {
        document.body.removeChild(brs[i]);
    }
    for (let i = 0; i < resultNodes.length; i++) {
        if (resultNodes[i] !== undefined) document.body.removeChild(resultNodes[i]);
    }
    resultNodes = [];
    height = +document.getElementById('matrix-height').value;
    width = +document.getElementById('matrix-width').value;
    nodes = [];
    brs = [];
    matrix = [];
    ref = [];
    pivotColumns = [];
    transposed = [];
    rref = [];

    if (!(height > 0 && width > 0)) return;

    for (let i = 0; i < height; i++) {
        nodes[i] = [];
        for (let j = 0; j < width; j++) {
            let node = document.createElement('input')
            node.type = 'number';
            node.id = `${i}-${j}`;
            node.className = 'element';
            nodes[i][j] = node;
            document.body.appendChild(node);
        }
        let br = document.createElement('br');
        brs[i] = br;
        document.body.appendChild(br);
    }

    //////////////////////////////////////////////////////////

}

function updateResults() {
    for (let i = 0; i < resultNodes.length; i++) {
        if (resultNodes[i] !== undefined) document.body.removeChild(resultNodes[i]);
    }
    resultNodes = [];

    for (let i = 0; i < nodes.length; i++) {
        matrix[i] = [];
        for (let j = 0; j < nodes[i].length; j++) {
            if (nodes[i][j].value === null || nodes[i][j] === "") return;
            matrix[i][j] = +nodes[i][j].value;
        }
    }
    let currentRow, currentColumn;
    [ref, currentRow, currentColumn] = REF(copy(matrix));
    [rref, pivotColumns] = RREF(copy(ref), currentRow, currentColumn);
    kernel = getCols(rref, invertArr(pivotColumns, width)).map(negativeRow);
    kernel = kernel.concat(identity(kernel[0].length));

    tableMatrix(fix(copy(ref)), "Row echelon form:", 0);
    tableMatrix(fix(copy(rref)), "Reduced row echelon form:", 1);
    tableMatrix(getCols(matrix, pivotColumns), "Col A:", 2);
    tableMatrix(kernel, "Kernel/null space: ", 3);
    tableMatrix(transpose(matrix), "Transposed:", 4);

}

function sleep(milliseconds) {
    const date = Date.now();
    let currentDate = null;
    do {
        currentDate = Date.now();
    } while (currentDate - date < milliseconds);
}

function tableMatrix(imatrix, label, index) {
    resultNodes[index] = document.createElement('div');
    resultNodes[index].innerText = label;
    let table = document.createElement('table');
    resultNodes[index].appendChild(table);
    for (let i = 0; i < imatrix.length; i++) {
        let tableRow = document.createElement('tr');
        for (let j = 0; j < imatrix[i].length; j++) {
            let element = imatrix[i][j];
            let value = document.createElement('td');
            value.innerText = +element === parseInt(element, 10) ? element : math.format(math.fraction(element));
            value.style.cssText = "width:40px;height:40px;text-align:center"
            tableRow.appendChild(value);
        }
        table.appendChild(tableRow);
    }
    document.body.appendChild(resultNodes[index]);
}

function negativeRow(row) {
    return row.map(function (x) {
        return -x
    });
}

function LCM(x, y) {
    if ((typeof x !== 'number') || (typeof y !== 'number'))
        return false;
    return (!x || !y) ? 0 : Math.abs((x * y) / GCD(x, y));
}

function GCD(x, y) {
    x = Math.abs(x);
    y = Math.abs(y);
    while (y) {
        let t = y;
        y = x % y;
        x = t;
    }
    return x;
}

function addRows(row1, row2) {
    for (let i = 0; i < row1.length; i++) {
        row1[i] += row2[i];
    }
    return row1;
}

function multRow(row, scalar) {
    return row.map(function (num) {
        return num * scalar
    });
}

function transpose(matrix) {
    let tmatrix = [];
    for (let i = 0; i < matrix.length; i++) {
        for (let j = 0; j < matrix[i].length; j++) {
            if (tmatrix[j] === undefined) tmatrix[j] = [];
            tmatrix[j][i] = matrix[i][j];
        }
    }
    return tmatrix
}

function fix(matrix) {
    for (let i = 0; i < matrix.length; i++) {
        for (let j = 0; j < matrix[i].length; j++) {
            matrix[i][j] = math.fraction(matrix[i][j]);
        }
    }
    return matrix;
}

function copy(x) {
    return JSON.parse(JSON.stringify(x));
}

function getCols(imatrix, cols) {
    let omatrix = [];
    for (let i = 0; i < imatrix.length; i++) {
        let row = [];
        for (let j = 0; j < cols.length; j++) {
            row.push(imatrix[i][cols[j]]);
        }
        omatrix[i] = row;
    }
    return omatrix;
}

function REF(matrix) {
    let currentColumn = 0;
    let currentRow = 1;
    while (currentColumn < (nodes[0].length - 1)) {
        if (currentRow >= nodes.length) {
            currentRow = currentColumn + 2;
            currentColumn++;
            continue;
        }
        if (matrix[currentRow][currentColumn] === 0) {
            currentRow++;
            continue;
        }
        let pivot = matrix[currentColumn][currentColumn];
        let zero = matrix[currentRow][currentColumn];
        let lcm = LCM(zero, pivot);
        let multPivotRow = multRow(matrix[currentColumn], (lcm / pivot));
        let multZeroRow = multRow(matrix[currentRow], (lcm / zero));
        matrix[currentRow] = addRows(multZeroRow, negativeRow(multPivotRow));
    }
    return [matrix, currentRow, currentColumn];
}

function RREF(imatrix, currentRow, currentColumn) {
    let omatrix = copy(imatrix);
    let pivotColumns = [];
    currentRow--;
    rref = copy(ref);
    for (let i = 0; i < rref.length; i++) {
        let pivot = rref[i].reduce(function (total, currentValue, currentIndex) {
            if (total === -1) {
                if (currentValue !== 0) return currentIndex;
            }
            return total;
        }, -1);
        if (pivot === -1) continue;

        pivotColumns.push(pivot);
    }
    for (let i = 0; i < pivotColumns.length; i++) {
        omatrix[pivotColumns[i]] = multRow(omatrix[pivotColumns[i]], 1 / omatrix[pivotColumns[i]][pivotColumns[i]]);
    }
    let pivotIndex = pivotColumns.length - 1;
    currentRow = pivotColumns[pivotIndex] - 1;
    currentColumn = currentRow - 1;
    do {
        console.log(pivotIndex);
        if (currentRow < 0) {

            pivotIndex--;
            currentColumn = pivotColumns[pivotIndex];
            currentRow = currentColumn - 1;
            continue;
        }
        if (omatrix[currentRow][pivotColumns[pivotIndex]] === 0) {
            currentRow--;
            continue;
        }
        omatrix[currentRow] = addRows(omatrix[currentRow], negativeRow(multRow(omatrix[pivotColumns[pivotIndex]], omatrix[currentRow][pivotColumns[pivotIndex]])));

    } while (currentRow >= 0 || pivotIndex !== 0)
    return [omatrix, pivotColumns];
}

function identity(n) {
    let imatrix = [];
    for (let i = 0; i < n; i++) {
        imatrix[i] = [];
        for (let j = 0; j < n; j++) {
            imatrix[i][j] = +(i === j);
        }
    }
    return imatrix;
}

function invertArr(arr, n) {
    let arr2 = [];
    for (let i = 0; i < n; i++) {
        if (!arr.includes(i)) arr2.push(i);
    }
    return arr2;
}