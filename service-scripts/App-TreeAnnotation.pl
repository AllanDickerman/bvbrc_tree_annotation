# The TreeAnnotation application.
use strict;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use File::Copy ('copy', 'move');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Bio::P3::Workspace::WorkspaceClientExt;
use Cwd;
use Phylo_Tree; # should be in lib directory

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;
my @default_genome_metadata_fields = (
        "genome_name", "family", "order", "isolation_country", "collection_year");
my @default_feature_metadata_fields = ("product", "accession");

our $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
if ($debug) {
    print STDERR "debug = $debug\n" if $debug;
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}
my $sequence_identifier_type; # feature_id or genome_id or user_specified

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#$data_url = "https://patricbrc.org/api" if $debug;
print STDERR "data_url=\n$data_url\n" if $debug;

my $wsClient = Bio::P3::Workspace::WorkspaceClientExt->new();
my $script = Bio::KBase::AppService::AppScript->new(\&annotate_tree, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";
    my $pf = {
	cpu => 1,
	memory => "1G",
	runtime => 60,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub annotate_tree {
    my ($app, $app_def, $raw_params, $params) = @_;

    print "Proc GeneTree build_tree ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my @outputs; # array of tuples of (filename, filetype)
    my $api = P3DataAPI->new;
    my $tmpdir = File::Temp->newdir( "/tmp/GeneTree_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
    my $original_wd = getcwd();
    chdir($tmpdir); # do all work in temporary directory
   
    my @feature_metadata_fields = @default_feature_metadata_fields;
    if (exists $params->{feature_metadata_fields}) {
        @feature_metadata_fields = @{$params->{feature_metadata_fields}};
    }
    #ensure that feature_id and genome_id are retrieved
    push @feature_metadata_fields, "feature_id" unless grep(/feature_id/, @feature_metadata_fields);
    push @feature_metadata_fields, "genome_id" unless grep(/genome_id/, @feature_metadata_fields);

    my @genome_metadata_fields = @default_genome_metadata_fields;
    if (exists $params->{genome_metadata_fields}) {
        @genome_metadata_fields = @{$params->{genome_metadata_fields}};
    }
    my $workspace_newick_file = $params->{newick_file};
    my $treeFile = "tree.nwk";
    $wsClient->download_file($workspace_newick_file, $treeFile);
    my $database_link = 'genome_id';

    print STDERR "tree_file $treeFile\n";
    if ($database_link) { # use system call to p3x-newick-to-phyloxml 
        my @command = ('p3x-newick-to-phyloxml', '-l', $database_link, '-g', join(',',@genome_metadata_fields), '-f', join(',', @feature_metadata_fields), $treeFile);
        print STDERR "execute system call: (as array):\n" . join(' ', @command), "\n";
        system(@command);
        my $phyloxml_file = $params->{output_file};
        push @outputs, [$phyloxml_file, "phyloxml"];
    }
    
    print STDERR '\@outputs = '. Dumper(\@outputs);
    my $output_folder = $app->result_folder();
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n" if $debug;
        if (0) { # for some reason this doesn't work
           $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
               (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
               $global_token);
        }
        else { # fall back to calling CLI
            my $ext = $1 if $ofile =~ /.*\.(\S+)$/;
            my @cmd = ("p3-cp", "-f", "-m", "${ext}=$type", $ofile, "ws:" . $app->result_folder);
            print STDERR "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        }
    }
    print STDERR "$tmpdir\n" if $debug;
    chdir($original_wd); # change back to the starting working directory
}

