
function initTable(table, taskid){
    table.DataTable({
        "columnDefs": [
            { "orderData": [ 0, 2 ], "targets": 0 },
            { "orderData": [ 0, 2 ], "targets": 1 },
            { "orderData": [ 2 ], "targets": 2 },
            { "orderData": [ 3 ], "targets": 3 },
            { "orderData": [ 4 ], "targets": 4 },
            { "orderData": [ 5 ], "targets": 5 },
            { "orderData": [ 6 ], "targets": 6 },
            { "orderData": [ 7 ], "targets": 7 },
            { "orderData": [ 8 ], "targets": 8 },
            { "orderData": [ 9 ], "targets": 9 },
            { "orderData": [ 10 ], "targets": 10 }
        ],
        "pageLength": 1000,
        // TODO handle multiple pages
        "initComplete": function(settings, json) {
            editTable(table);
            addAlignments(table, taskid)
        }
    });
}

function showTooltips() {
    $('body').tooltip({
        selector: '.alignment .residue',
        delay: { "show": 200, "hide": 50 },
        html: true,
        title: function(){
            var labels = []
            var ref = $(this).data('ref')
            var mut = $(this).text()
            if (ref == '-') {
                labels.push('insertion')
            } else if (mut == '-') {
                labels.push('deletion')
            }
            var pos = $(this).data('display-pos')
            labels.push('position: '+ pos)
            labels.push('reference: <b>'+$(this).data('ref')+'</b>')
            labels.push('traces: '+$(this).data('options'))
            labels.push("Click to show traces")
            return labels.join('<br>')
        }
    });
}

// Add alignment as child to corresponding results row
function addAlignments(table, taskid) {
    table.find('tr').each(function () {
        var alignment = $(this).find('.alignment').removeClass('collapse')
        var row = table.DataTable().row($(this));
        row.child(alignment);
        // underlining codons
        alignment.find('.residue').mouseover(codonMouseover);
        alignment.mouseout(codonMouseout);
        // display trace viewer on click
        alignment.find('.residue').click(showTraceViewer);


        function codonMouseover(){
            alignment.find('.codon-ref').removeClass('codon-ref')
            alignment.find('.codon-new').removeClass('codon-new')

            var res = $(this);
            var num_ref = res.data('codon-num-ref')
            var ref_base = res.data('ref')
            var num_new = res.data('codon-num-new')
            var new_base = res.text()
            var featureID = res.data('feat')

            if (ref_base != '-' && typeof num_ref != "undefined") {
                res.addClass('codon-ref')
                res.nextUntil().filter(function () {
                        return num_ref == $(this).data('codon-num-ref') && featureID == $(this).data('feat');
                    }).addClass('codon-ref')
                res.prevUntil().filter(function () {
                        return num_ref == $(this).data('codon-num-ref') && featureID == $(this).data('feat');
                    }).addClass('codon-ref')
            }

            if (new_base != '-' && typeof num_new != "undefined") {
                res.addClass('codon-new')
                res.nextUntil().filter(function () {
                        return num_new == $(this).data('codon-num-new') && featureID == $(this).data('feat');
                    }).addClass('codon-new')
                res.prevUntil().filter(function () {
                        return num_new == $(this).data('codon-num-new') && featureID == $(this).data('feat');
                    }).addClass('codon-new')
            }
        }

        function codonMouseout(){
            alignment.find('.codon-ref').removeClass('codon-ref')
            alignment.find('.codon-new').removeClass('codon-new')
        }

    });
    table.on('click', '.ref-button', showAlignment)

    // Display or hide alignment row on reference ID button click
    function showAlignment() {
        $(this).closest('tr').each(function(){
            var row = table.DataTable().row( $(this) );
            if ( row.child.isShown() ) {
                row.child.hide();
            } else {
                row.child.show();
            }
        })
    };

    function showTraceViewer(){
        var res = $(this);
        var pos = res.data('position')
        var alignment = res.closest('.alignment')
        var alignmentIndex = alignment.data('alignment')

        if (pos == null) {
            return;
        }

        var errorElement = $('#traceError')
        var titleElement = $('.modal-title')
        titleElement.text('Loading...')
        $("#traceModal svg").remove();
        $('#traceModal').modal('show');
        urlstr = "/trace/".concat(taskid.toString()).concat("/");

        $.ajax({
            type: "GET",
            url: urlstr.concat(alignmentIndex.toString()),
            success: function(response) {
                errorElement.text('')
                titleElement.html("Traces for Reference ".concat(response.refid))
                var viewer = new TraceViewer(response);
                viewer.show(pos);
            },
            error: function(xhr, ajaxOptions, thrownError) {
                titleElement.text('Error')
                errorElement.html(xhr.responseText)
            }
        });
    }
}


function tableColors() {
    $('td[data-color-percentage]').each(function(){
        var percentage = $(this).data('color-percentage')
        var bin = $(this).data('color-bin') + 0
        var bins = [
            'rgb(240,  70, 70)',
            'rgb(255,  83, 83)',
            'rgb(255, 128, 57)',
            'rgb(255, 172, 51)',
            'rgb(250, 193, 81)',
            'rgb(150, 198,106)',
            'rgb(95,  200, 95)'
        ]
        $(this).css('background', 'linear-gradient(90deg, '+bins[bin]+' '+percentage+'%, transparent '+percentage+'%)')
        $(this).css('border-right', '1px solid #cccccc')
    });
}

function editTable(table){
    var collapseButton = $('<button type="button" class="btn btn-primary" value="Expand all" id="collapseButton">Expand all</button>').click(function() {
        table.find('tr').each(function(){
            var row = table.DataTable().row( $(this) );
            if ( $('#collapseButton').val() == "Collapse all" ) {
                row.child.hide();
            } else {
                row.child.show();
            }
        })
        if ($(this).text() == "Collapse all"){
            $(this).text("Expand all")
            $(this).val("Expand all")
        } else {
            $(this).text("Collapse all")
            $(this).val("Collapse all")
        }
    });
    $('.dataTables_length').prepend(collapseButton);

    // Prettify search field
    $('.dataTables_filter').addClass('form-inline').find('input[type=search]').addClass('form-control form-control-sm');

}


