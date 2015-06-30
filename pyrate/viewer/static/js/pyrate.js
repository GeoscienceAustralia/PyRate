var currentLayer = undefined;
// assumes the following are set (which currently is done in index.html by the
// template engine.
// - MISSING_VALUE
// - PYRATE_BASE_URL



// @brief Get the timeseries for a location and plot it.
//
// @param coordinate The coordinates (in the projection used by the map) of the
// location clicked.
//
// @param dataPlotter An instance of a data plotter (as returned by the function
// *Plotter* below).
function getTimeSeries(coordinate, dataPlotter) {
    var success = function(resp) {
        if(resp.exception != undefined) {
            alert(resp.exception);
            return;
        }
        dataPlotter.showTimeSeries(resp);
    }
    $.ajax({
        url: PYRATE_BASE_URL + '/ts',
        data: {x: coordinate[0], y:coordinate[1]},
        error: function(resp) { alert('Failed to get time series: ' + resp.toString()); },
        success: success,
        context: dataPlotter,
        type: 'GET',
        async: true
    });
}



// @brief Setup the the dropdown list used to select the file to show.
//
// This attaches a function to the dropdown list which calls *callback* when
// the selected item in the list changes. It also makes a call to the server
// to get the list of files to populate the dropdown list with.
//
// @param dropdownListSelector JQuery selector for the dropdown list.
//
// @param callback function to be called when the selection changes. This must
// take a single argument which is the name of the file selected.
function setupFileSelect(dropdownListSelector, callback) {
    var elem = $(dropdownListSelector);
    elem.off('change');

    // this function is used to call *callback* when the selected item in the
    // dropdown list is changed.
    var imageFetcher = function() {
        var fileName = elem.find(':selected').text();
        if(fileName != '') {
            callback(fileName);
        }
    }

    // attach the above function to the change event of the dropdown list
    elem.change(imageFetcher);

    // function to populate the dropdown list upon success of the ajax call
    var success = function(resp) {
        if(resp.exception != undefined) {
            alert(resp.exception);
            return;
        }
        var htmlStr = '<option></option>';
        $.each(resp, function(index, fileName) {
            htmlStr += '<option>' + fileName + '</option>';
        });
        elem.html(htmlStr);
        elem.prop('disabled', false);
    };

    // make the ajax call to the server.
    $.ajax({
        url: PYRATE_BASE_URL + '/files',
        error: function(resp) { alert('Failed to get list of files.'); },
        success: success,
        type: 'GET',
        async: true // not sure if I want this or not... just playing for now.
    });

    return imageFetcher;
}



// @brief Get the image for a tif file.
//
// This makes an ajax call back to ther server to retrieve the image for
// *fileName*. Upon success of that call, a new OpenLayers layer for the
// retrieved image is created and passed to *callback*.
//
// @param fileName The name of the tif file to get an image for.
//
// @param callback A callback which will be called with the new layer as the
// first argument, and the response from server as a second.
//
// @param discreteColors A boolean determining whether to apply a discrete
// or continuous color scheme.
function getImageForFileName(fileName, callback, discreteColors) {
    var success = function(result) {
        if(result.exception != undefined) {
            alert(result.exception);
            return;
        }
        // add the layer to the map
        var newLayer = new ol.layer.Image({
            source: new ol.source.ImageStatic({
                url: result.url,
                projection: result.projection,
                imageExtent: result.imageExtent,
                imageSize: result.imageSize
            })
        });
        callback(newLayer, result);
    };
    $.ajax({
        url: PYRATE_BASE_URL + '/image/' + fileName,
        error: function(resp) { alert('Failed to retrieve image for ' + fileName + '.'); },
        data: { discrete : discreteColors },
        success: success,
        type: 'GET',
        async: false
    });
}



// @brief Make a function that can be used to retrieve an image for a file name.
//
// The function returned is called with a single tif file name and populates
// the map and legend upon sucess.
//
// @param map The Openlayers map instance.
//
// @param discreteChecker A no argument function that when called will return
// whether to request a discrete or continuous color scheme for the returned
// layer.
//
// @param topLayer A layer that will be placed on top of all other others after
// the result image placed on the map. This is used to keep the layer on which
// transects are drawn on top.
//
// @colorTableDivSelector The JQuery selector for the div which holds the legend
// for the map.
function makeRegionImageFetcher(map, discreteChecker, topLayer, colorTableDivSelector) {

    // function to populate the color table
    var generateColorTable = function(colorTable) {
        var newReplaced = $('<div></div>');
        newReplaced.attr('class', 'color-entry');
        for(var level=-1; level!=colorTable.length; ++level) {
            var div = $('<div>');
            var spn = $('<span></span>');
            if(level == -1) {
                spn.append('Legend');
            } else {
                var vals = colorTable[level][1];
                var rgb = 'rgb('+vals[0]+','+vals[1]+','+vals[2]+')';
                div.css('background-color', rgb);
                if(discreteChecker()) {
                    if(level == 0) {
                        spn.append('< ' + parseFloat(colorTable[level]).toFixed(2));
                    } else if(level == (colorTable.length - 1)) {
                        spn.append('> ' + parseFloat(colorTable[level]).toFixed(2));
                    } else {
                        spn.append(
                            parseFloat(colorTable[level-1]).toFixed(2) +
                            ' ~ ' +
                            parseFloat(colorTable[level]).toFixed(2)
                        );
                    }
                } else {
                    spn.append(parseFloat(colorTable[level]).toFixed(2));
                }
            }
            div.append(spn);
            newReplaced.append(div);
        }
        $(colorTableDivSelector).html(newReplaced);
    };

    // function to put the result image on the map and return *topLayer* to the
    // top.
    return function(fileName) {
        var callback = function(newLayer, result) {
            map.removeLayer(currentLayer);
            currentLayer = newLayer;
            map.addLayer(currentLayer);
            currentLayer.setVisible(true);
            map.getLayers().setAt(map.getLayers().getArray().length, topLayer);
            // make the color table
            if(result.colorTable != undefined) {
                $(colorTableDivSelector).show();
                generateColorTable(result.colorTable);
            } else {
                $(colorTableDivSelector).hide();
            }
        };
        getImageForFileName(fileName, callback, discreteChecker());
    };
}



// @brief Create a drawing layer for the map.
//
// Most of this was adapted from the example code at
// http://openlayers.org/en/v3.6.0/examples/line-arrows.html.
//
// @param map An openlayers map instance.
//
// @param dataPlotter An object used for plotting transects. It must have a
// function attribute *showTransect* which will be passed the response of an
// ajax call to */transect* (see the Python function *pyrate.viewer.web.transect*).
// At present, *dataPlotter* is created by a call to the function *Plotter* below.
//
// @param epsilonTextfieldSelector JQuery selector for the text field
// containing the value of epsilon (which is used for determining the weights
// to use for neighbouring cells - see the Python function
// *pyrate.viewer.web.smoother* for how it is used).
//
// @param maxNeighborhoodTextfieldSelector JQuery selector for the text field
// containing the value specifying the maximum width of the neighbourhood. See
// the Python function *pyrate.viewer.web.smoother* for how it is used.
//
// Returns a object with 3 attributes: *drawer*, *layer* and *update*. *drawer*
// is an instance of *ol.interaction.Draw* which is used to draw transects.
// When the end point of the transect is drawn, an ajax call is made to get
// the data for the transect. *layer* is the layer on which the transects are
// drawn and hence should be kept on top of all other layers. *update* is a
// function that can be called to update the parameters used in the calculation
// of the trasects.
function createDrawingLayer(map, dataPlotter, epsilonTextfieldSelector, maxNeighborhoodTextfieldSelector) {

    // holds the features for a vector layer.
    var source = new ol.source.Vector();

    // copied from http://openlayers.org/en/v3.6.0/examples/line-arrows.html
    var styleFunction = function(feature, resolution) {
        var geometry = feature.getGeometry();
        var styles = [
            // linestring
            new ol.style.Style({
                stroke: new ol.style.Stroke({
                    color: '#ffcc33',
                    width: 2
                })
            })
        ];

        geometry.forEachSegment(function(start, end) {
            var dx = end[0] - start[0];
            var dy = end[1] - start[1];
            var rotation = Math.atan2(dy, dx);
            // arrows
            styles.push(new ol.style.Style({
                geometry: new ol.geom.Point(end),
                image: new ol.style.Icon({
                    src: '/static/arrow.png',
                    anchor: [0.75, 0.5],
                    rotateWithView: false,
                    rotation: -rotation
                })
            }));
        });

        return styles;
    };

    // the vector layer transects are drawn on
    var transectLayer = new ol.layer.Vector({
        source: source,
        style: styleFunction
    });

    // the tool that draws the transects
    var drawer = new ol.interaction.Draw({
        source: source,
        type: /** @type {ol.geom.GeometryType} */ ('LineString'),
        maxPoints: 2
    });

    // function that extracts the relevant parameters from the page and makes
    // the ajax request to the server.
    var drawTransects = function() {
        var epsilon = $(epsilonTextfieldSelector).val();
        var maxNeighborhood = $(maxNeighborhoodTextfieldSelector).val();
        $.ajax({
            url: PYRATE_BASE_URL + '/transect',
            error: function(resp) { alert('Failed to retrieve transect'); },
            data: {
                geom:JSON.stringify({start:start, end:end}),
                epsilon:epsilon,
                maxNeighborhood:maxNeighborhood
            },
            success: dataPlotter.showTransect,
            context: dataPlotter,
            type: 'GET',
            async: true
        });
    };

    // the start and end of the transect. kept here so that the transect can be
    // redrawn when, for example, the smoothing parameters are changed.
    var start = null;
    var end = null;

    // wire up the action to be performed when the second (final) point in a
    // transect is drawn.
    drawer.on('drawend', function(evt) {
        self = this;

        // the feature that has just been drawn
        var feature = evt.feature;

        // remove all features except the one just drawn
        self.getSource().forEachFeature(
            function(f) { if(f != feature) this.removeFeature(f); },
            self.getSource()
        );

        // the the start and end points of the transect
        var geom = feature.getGeometry();
        start = geom.getCoordinates()[0];
        end = geom.getCoordinates()[1];

        drawTransects();
    }, transectLayer);

    return {
        drawer: drawer,
        layer: transectLayer,
        update: function() {
            if(start && end) {
                drawTransects();
            } else {
                alert("please draw a transect");
            }
        }
    };
}



// @brief Make an object responsible for managing the plots.
//
// @param mainPlotDivSelector JQuery selector for the div in which to place
// the 'main' plot.
//
// @param weightsDivSelector JQuery selector for the div in which to place
// weights plot.
//
// @param legendDivSelector JQuery selector for the div in which to place the
// legend for the main plot.
function Plotter(mainPlotDivSelector, weightsDivSelector, legendDivSelector) {

    // turn arrays *x* and *y* into a suitable data structure for plotting in flot.
    function makeSeries(datatype, x, y) {
        if(x.length != y.length) {
            return false;
        }
        var xyData = [];
        for(var i=0; i<x.length; ++i) {
            var value = y[i];
            if(value != MISSING_VALUE) { // assumes that integer like comparision will work.
                xyData.push([x[i], y[i]]);
            }
        }
        return {
            datatype: datatype,
            data: xyData
        }
    }

    // make the label for a dataset (used in legends)
    function makeDatasetLabel(data) {
        return data.datatype;
    }

    // make a series suitable for plotting in flot.
    function makeDataset(dataset, index, color) {
        var tag = makeDatasetLabel(dataset);
        return {
            index: index,
            tag: tag,
            color: color,
            label: tag,
            data: dataset.data,
            shadowSize: 0.0
        };
    }

    // parse a date from a string like *yyyy-mm-yy*.
    function parseDate(input) {
        var parts = input.split('-');
        return new Date(parts[0], parts[1]-1, parts[2]); // Note: months are 0-based
    }

    dataPlotter = {
        elem: $(mainPlotDivSelector),
        weights_elem: $(weightsDivSelector),
        the_plot: null,
        the_weights_plot: null,
        hoverPrevious: null,

        // create a popup for a point in a data series in the main plot.
        drawToolTip: function(x, y, content) {
            $("<div id='chart-tooltip'>" + content + "</div>").css({
                'z-index': 1001,
                position: "absolute",
                display: "none",
                top: y + 5,
                left: x + 5,
                border: "1px solid #fdd",
                padding: "2px",
                color: "white",
                'background-color': "#222",
                opacity: 0.80,
                'font-size': 'small',
                'border-radius': '5px'
            }).appendTo("body").fadeIn(200);
        },

        // make a callback for when a data point in the main plot is clicked
        // results in a call to drawToolTip
        makePlotHoverCallback: function(self) { return function(evt, pos, item) {
            if (item) {
                if (self.hoverPrevious != item.dataIndex) {
                    self.hoverPrevious = item.dataIndex;
                    $("#chart-tooltip").remove();
                    var d = new Date(item.datapoint[0]).toISOString().slice(0,11),
                        y = item.datapoint[1].toFixed(3);
                    var label = "Date: " + d + "<br/>Value:" + y + "<br/>";
                    self.drawToolTip(item.pageX, item.pageY, label);
                }
            } else {
                $("#chart-tooltip").remove();
                self.hoverPrevious = null;
            }
        }},

        // flot options for the main plot
        options:{
            series: {
                lines: {
                    show: true,
                    lineWidth: 0.5
                },
                points: {
                    show: true,
                    radius: 1
                }
            },
            grid: {
                hoverable: true,
                clickable: true,
                show: true
            },
            yaxis: {
                show: true,
                zoomRange: false,
                panRange: false
            },
            zoom:  {
                interactive: true
            },
            pan: {
                interactive: true
            },
            legend: {
                show: true,
                container: $(legendDivSelector)
            }
        },

        // flot options for the weights plot
        weights_options:{
            series: {
                bars: {
                    show: true
                }
            },
            bars: {
                align: "center",
                barWidth: 0.5
            },
            grid: {
                hoverable: true,
                clickable: true,
                show: true
            },
            xaxis: {
                show: true,
            },
            yaxis: {
                show: true,
                zoomRange: false,
                panRange: false
            },
            zoom:  {
                interactive: false
            },
            pan: {
                interactive: false
            }
        },

        // hide one or both plots
        hide: function() {
            if(this.the_plot) {
                this.the_plot.getPlaceholder().hide();
                if(this.the_weights_plot) {
                    this.the_weights_plot.getPlaceholder().hide();
                }
            }
        },

        // show one or both plots
        show: function() {
            if(this.the_plot) {
                this.the_plot.getPlaceholder().show();
            }
            if(this.the_weights_plot) {
                this.the_weights_plot.getPlaceholder().show();
            }
        },

        // redraw the plots
        reDraw: function() {
            if(this.the_plot) {
                this.show();
                this.the_plot.resize();
                this.the_plot.setupGrid();
                this.the_plot.draw();
                if(this.the_weights_plot) {
                    this.the_weights_plot.resize();
                    this.the_weights_plot.setupGrid();
                    this.the_weights_plot.draw();
                }
            }
        },

        // setup and draw the plots
        setupAndDraw: function() {
            $("#chart-tooltip").remove();
            if(this.datasets.length === 0) {
                this.hide();
            } else {
                this.show();
                if(this.isTransect) {
                    this.options.xaxis = { show: true };
                    this.the_plot = $.plot(this.elem, this.datasets, this.options);
                    this.the_weights_plot = $.plot(this.weights_elem, this.kernel_weights, this.weights_options);
                } else {
                    this.options.xaxis = {show: true, mode: "time", timeformat: "%Y/%m/%d", rotateTicks: 45};
                    this.the_plot = $.plot(this.elem, this.datasets, this.options);
                    this.the_weights_plot = null;
                }
                this.the_plot.getPlaceholder().bind("plotclick", this.makePlotHoverCallback(this));
                $('<div class="chart-button" style="right:5px;top:50px">zoom<br/>in</div>')
                    .appendTo(this.elem).click(function (e) {
                        e.preventDefault();
                        dataPlotter.the_plot.zoom();
                });
                $('<div class="chart-button" style="right:5px;top:100px">zoom<br/>out</div>')
                    .appendTo(this.elem).click(function (e) {
                        e.preventDefault();
                        dataPlotter.the_plot.zoomOut();
                });
                $('<div class="chart-button" style="right:5px;top:150px">zoom<br/>reset</div>')
                    .appendTo(this.elem).click(function (e) {
                        e.preventDefault();
                        // there must be a nicer way to do this!
                        dataPlotter.the_plot = null;
                        dataPlotter.setupAndDraw();
                });
            }
        },

        // Show a time series. This is intended to be used as a calback for an
        // ajax request for a time series.
        showTimeSeries: function(resp) {
            this.datasets = [];
            for(var i=0; i<resp.times.length; ++i) {
                resp.times[i] = parseDate(resp.times[i]);
            }
            this.datasets.push(makeDataset(makeSeries(
                'data',
                resp.times,
                resp.data
            ), "#3399FF"));
            this.isTransect = false;
            this.setupAndDraw();
        },

        // show the data for a transect. This is intended to be used as a
        // calback for an ajax request for a time series.
        showTransect: function(resp) {
            this.datasets = [];
            this.kernel_weights = [];
            this.kernel_weights.push(makeDataset(makeSeries('weights', resp.weights.x, resp.weights.y), "#3399FF"));
            for(var i=0; i<resp.transects.length; ++i) {
                var y = resp.transects[i]['data'];
                var x = [];
                var color = $.Color(33, 99, i*255/resp.length).toHexString(false);
                for(var j=0; j<y.length; ++j) x.push(j);
                this.datasets.push(makeDataset(makeSeries(resp.transects[i]['date'], x, y), color));
            }
            this.isTransect = true;
            this.setupAndDraw();
        }
    };

    return dataPlotter;
}
