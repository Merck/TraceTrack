function TraceViewer(data) {

    /**
    DATA
    **/

    // list of position dictionaries with:
    //      cov: coverage (int)
    //      mut: type of mutation (str)
    //      ref: corresponding position in unaligned reference (int)
    this.positions = data.positions

    // Extracted list where alignedRefPositions[mutatedPos] = refPos
    this.alignedRefPositions = this.positions.map(d => d.ref)

    // list of trace sequence dictionaries with:
    //      id: ID based on trace filename
    //      alignedPositions: list where alignedPositions[mutatedPos] = traceSeqPos
    //      sequence: original unaligned trace sequence (str)
    //      locations: list of trace x-coordinates of each sequence base
    //      traces: list of trace value dictionaries with:
    //          base: A, C, T or G
    //          values: list of trace values
    this.sequences = data.sequences

    // mutated result sequence (str)
    this.mutated = data.mutated

    // original unaligned reference sequence (str)
    this.reference = data.reference

    // translation of reference sequence (str)
    this.translation = data.translation

    // ID of reference sequence
    this.refid = data.refid

    /**
    PLOT PROPERTIES
    **/

    this.refRange = 10
    this.traceRange = 112
    this.margin = {top: 20, right: 20, bottom: 40, left: 40};
    this.totalWidth = +d3.select('#tracePlot').style('width').slice(0, -2);
    if (!this.totalWidth) {
        // TODO

        this.totalWidth = 780
    }
    this.width = this.totalWidth - this.margin.left - this.margin.right;
    this.height = 75;
    var viewer = this

    /**
     DRAW NAVIGATION BAR
    **/

    this.updateNavArrow = function(mutatedPos) {  // Add interactivity
        var mutBase = viewer.mutated[mutatedPos]
        var refPos = viewer.positions[mutatedPos].ref
        var textPos = viewer.positions[mutatedPos].pos
        var text = ''
        if (mutatedPos < viewer.positions.length - 1 && refPos == viewer.positions[mutatedPos + 1].ref) {
            text += textPos + ' insert: '+mutBase
        } else {
            var refBase = viewer.reference[refPos]
            text += (refPos + 1) + " " + refBase
            if (mutBase != refBase) {
                if (mutBase == '-') {
                    text += ' (deletion)'
                } else {
                    text += ' (mut '+mutBase+')'
                }
            }
        }

        viewer.navArrow
            .attr("transform", "translate(" + viewer.xScaleNav(mutatedPos) + ",10)")
            .select("text")
            .text(text)
    }

    this.xScaleNav = d3.scaleLinear()
        .domain([0, this.positions.length])
        .range([15, this.totalWidth-15]); // output

    // add navigator svg
    this.navPlot = d3.select('#navPlot').append("svg")
        .attr("width", this.totalWidth)
        .attr("height", 80)
      .append("g")
        .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")")

    this.navPlot.on("mouseout", function(d, i) {
        viewer.updateNavArrow(viewer.shownPos)
    })

    this.navPlot.append("text")
        .text("Navigator")
        .attr("x", this.width/2).attr("y", -5)
        .attr("class", "plotTitle");

    // add navigator rectangles
    var navRectWidth = Math.max(1, this.totalWidth / this.positions.length + 0.1)
    this.navPlot
        .append("g")
        .attr("transform", "translate(0, 15)")
        .selectAll("rect")
        .data(this.positions)
      .enter()
        .append("rect")
        .attr("width", navRectWidth)
        .attr("height", function(d, i){ return viewer.mutated[i] == viewer.mutated[i].toUpperCase() ? 30 : 24 })
        .attr("y", function(d, i){ return viewer.mutated[i] == viewer.mutated[i].toUpperCase() ? 0 : 3 })
        .attr("x", function(d, i){ return viewer.xScaleNav(i) - navRectWidth/2 })
        .attr("class", function(d){ return d.mut + ' cov' + Math.min(d.cov, 4) })
        .on("mouseover", function(d, i) { viewer.updateNavArrow(i) })
        .on("click", function(d, i) {
            viewer.show(i)
        })

    this.navArrow = viewer.navPlot.append("g")
        .attr("class", "navArrow")

    this.navArrow.append("text")

    this.navArrow.append("line")
        .attr("x1", 0)
        .attr("y1", 2)
        .attr("x2", 0)
        .attr("y2", 35)
        .style("stroke-width", 1)
        .style("stroke", "black")
        .style("fill", "none")


    /**
     DRAW REFERENCE
    **/

    this.xScaleRef = d3.scaleLinear()
        .range([0, this.width]); // output

    // Add SVG to the page to "ref" div
    this.refPlot = d3.select("#refPlot").append("svg")
        .attr("width", this.totalWidth)
        .attr("height", 75)
      .append("g")
        .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")")

    // add text label with file name
    this.refPlot.append("text")
        .attr("x", this.width/2).attr("y", -10)
        .text('Reference '+this.refid)
        .attr("class", "plotTitle");

    var refTicks = this.refPlot.append('g')
        .selectAll('g.refTick')
        .data(this.reference.split(''))
      .enter()
        .append('g')
        .attr('class', function(d, i){ return 'refTick tick tick_'+d.toUpperCase() });

    refTicks.insert('rect')
        .attr('transform', 'translate(-10 -10)')
        .attr('width', 20)
        .attr('height', 20)
        .attr('rx', 4)
        .on("click", function(d, refPos){
            var mutatedPos = viewer.alignedRefPositions.lastIndexOf(refPos)
            viewer.show(mutatedPos)
         })

    refTicks.insert('text')
        .attr("x", 0).attr("y", -20)
        .attr("dominant-baseline", "central")
        .text(function(d, i){ return i + 1 })
        .attr("class", "refPosition");

    refTicks.insert('text')
        .attr("x", 0).attr("y", -1)
        .attr("dominant-baseline", "central")
        .text(function(d, i){ return d });


    amino_acid = function(d, i) {
        return viewer.translation[i]
    }

    refTicks.insert('text')
        .attr("x", 0).attr("y", 20)
        .attr("dominant-baseline", "central")
        .text(amino_acid)
        .attr("class", "refTranslation");


    // Add SVG to the page to "ref" div
    this.navArrowL = d3.select("#refPlot").append("svg")
        .attr("width", 20)
        .attr("height", 20)
        .style("z-index", 1000)
        .style("position", "relative")
        .attr("transform", "translate(" + -5 + "," + 0 + ")")

    this.navArrowL.append('text')
        .text("◀")
        .attr("x", 10).attr("y", 15)
        .attr("class", "plotTitle")
        .style("cursor", "pointer")
        .on("click", function(){
            viewer.show(viewer.shownPos - (viewer.refRange*2))
         })

     this.navArrowR = d3.select("#refPlot").append("svg")
        .attr("width", 20)
        .attr("height", 20)
        .attr("transform", "translate(" + 760 + "," + 0 + ")")
        .style("z-index", 1000)
        .style("position", "relative")

     this.navArrowR.append('text')
        .text("▶")
        .attr("x", 10).attr("y", 15)
        .attr("class", "plotTitle")
        .style("cursor", "pointer")
        .on("click", function(){
            viewer.show(viewer.shownPos + (viewer.refRange*2))
         })

    /**
     DRAW TRACE PLOTS
    **/

    this.yScale = d3.scaleLinear()
        .domain([0, getMaxTrace(this.sequences)]) // input
        .range([this.height, 0]); // output

    // append svg for each trace sequence
    this.plots = d3.select('#tracePlot').selectAll("svg")
        .data(this.sequences)
      .enter()
        // enter an svg element for each tracefile
        .append("svg")
        .attr("width", this.totalWidth)
        .attr("height", this.height + this.margin.top + this.margin.bottom)
        .append("g")
        .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")")
        .attr("class", "main")

    // initialize each trace plot
    this.plots.each(function(d, i){
        // i is the trace file number (index in the list of dictionaries) d is dict
        var svg = d3.select(this);

        // add clipping rectangle mask
        svg.append("clipPath")
          .attr("id", "clipid")
        .append("rect")
          .attr("x", 0)
          .attr("y", -5)
          .attr("width", viewer.width)
          .attr("height", viewer.height+30);

        // add empty trace paths
        var paths = svg.append("g")
            .attr("class", "paths")
            .selectAll('path')
            .data(d.traces)
          .enter()
            .append('path')
            .attr('class', function(d){ return 'tracePath line_'+d.base })
            .attr("clip-path", "url(#clipid)");

        // add dashed line in center
        svg.append("line")
            .attr("x1", viewer.width / 2)
            .attr("y1", 0)
            .attr("x2", viewer.width / 2)
            .attr("y2", viewer.height)
            .style("stroke-dasharray", ("3, 3"))
            .style("stroke-width", 1)
            .style("stroke", "#aaaaaa")
            .style("fill", "none");

        // append empty x axis group
        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + viewer.height + ")");

        // append empty y axis group
        svg.append("g")
            .attr("class", "y axis")
            .call(d3.axisLeft(viewer.yScale));

        // add text label with file name
        svg.append("text")
            .attr("x", viewer.width/2).attr("y", 0)
            .text('Trace '+d.id)
            .attr("class", "plotTitle");

    });

    var transition = d3.transition()
        .delay(0)
        .duration(400)

    this.show = function(mutatedPos){

        this.shownPos = mutatedPos
        this.updateNavArrow(mutatedPos)

        var refPos = this.positions[mutatedPos].ref
        var opacity = 1
        if (mutatedPos < this.positions.length - 1) {
            // for insertions, center between two residues
            var nextRefPos = this.positions[mutatedPos+1].ref
            if (refPos == nextRefPos) {
                refPos = refPos - 0.5
                opacity = 0.5
            }
        }

        this.xScaleRef
            .domain([ refPos - this.refRange, refPos + this.refRange ])

        this.refPlot.selectAll("g.refTick")
            .transition(transition)
            .style("opacity", opacity)
            .attr("transform", function(d, i) { return "translate(" + viewer.xScaleRef(i) + ",25)" })


        this.plots.each(function(d, i){       // i is the trace file number (index in the list of dictionaries) d is dict
            // for each trace file edit svg and bind data
            var svg = d3.select(this)
            var tracePos = d.alignedPositions[mutatedPos]
            var nextTracePos = d.alignedPositions[mutatedPos+1]
            var traceCenter = d.locations[tracePos]
            var opacity = 1

            if (tracePos == nextTracePos) {
                // gap in trace, show with half opacity and in middle of previous two positions
                if (tracePos == 0) {
                    traceCenter = -10
                } else if (tracePos == d.locations.length) {
                    traceCenter = d.locations[d.locations.length-1] + 10
                } else {
                    traceCenter = (d.locations[tracePos] + d.locations[tracePos-1]) / 2
                }
                opacity = 0.5
            }

            var traceLength = d.traces[0].values.length
            var traceDrawMargin = 20
            var traceStart = traceCenter - viewer.traceRange
            var traceSliceStart = Math.max(0, traceStart - traceDrawMargin)
            var traceEnd = traceCenter + viewer.traceRange

            var xScale = d3.scaleLinear()
                .domain([ traceStart, traceEnd ]) // input
                .range([0, viewer.width]) // output

            // filter out base positions between start and end -> axis tick labels
            var ticks = d.locations.filter(function(number) {
              return number > traceStart && number <= traceEnd;
            });

            var seqStart = d.locations.indexOf(ticks[0]);
            var seqEnd = seqStart + ticks.length
            var tickLabels = d.sequence.split("").slice(seqStart, seqEnd)

            var xAxisGenerator = d3.axisBottom(xScale)
                .tickValues(ticks)
                .tickFormat((d,i) => tickLabels[i]);

            var line = d3.line()
                .x(function(d, i) { return xScale(traceSliceStart + i); }) // set the x values for the line generator
                .y(function(d) { return viewer.yScale(d); }) // set the y values for the line generator
                .curve(d3.curveMonotoneX) // apply smoothing to the line

            // trace line plot
            var paths = svg.selectAll('path.tracePath')
                .style("opacity", opacity)
                .attr("d", function(d){
                    return line(d.values.slice(traceSliceStart, traceEnd + traceDrawMargin))
                 });

            // Call the x axis in a group tag
            var xAxis = svg.select("g.x.axis")
                .call(xAxisGenerator)

            xAxis.selectAll('g.tick')
                .style("opacity", opacity)
                .each(function(unused, i){
                    // note: This only works if there is no transition
                    // because extra elements at the edge need to be added and animated
                    var g = d3.select(this);
                    var base = g.select('text').text()
                    var tick_class = "tick_uni"
                    if (base.toUpperCase() == "N"){
                        tick_class = "tick_N"
                    }

                    if (!g.select('rect').size()) {
                        var mutatedPos = d.alignedPositions.lastIndexOf(seqStart + i)
                        g.insert('rect', ':first-child')
                            .attr('transform', 'translate(-8 5)')
                            .attr('width', 16)
                            .attr('height', 17)
                            .attr('rx', 4)
                            .on("click", function(){
                                viewer.show(mutatedPos)
                             })
                    var refPos = viewer.alignedRefPositions[mutatedPos]
                    if (refPos && base.toUpperCase() != viewer.reference[refPos].toUpperCase()) {
                        tick_class = 'tick_'+base.toUpperCase()
                    }
                    g.attr('class', 'traceTick tick '+tick_class)

                    var prevMutPos = d.alignedPositions.lastIndexOf(seqStart + i - 1)
                    var prevRefPos = viewer.alignedRefPositions[prevMutPos]

                    if (refPos > prevRefPos + 1) {
                        g.insert('line', ':first-child')
                            .attr('transform', 'translate(-26 9)')
                            .attr('x1', 1)
                            .attr('y1', 6)
                            .attr('x2', 10)
                            .attr('y2', 6)
                            .attr('stroke-width', 2)
                            .attr('stroke', 'red')
                    }

                    }
                })

        });
    }
}

function getMaxTrace(sequences){
    var ys = []
    for (dict of sequences) {
        for (var i=0; i<4; i+=1) {
            ys.push(d3.max(dict.traces[i].values));
        }
    }
    return d3.max(ys);
}