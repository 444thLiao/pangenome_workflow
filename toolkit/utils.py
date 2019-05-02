from subprocess import check_call
def run_cmd(cmd, dryrun=False, **kwargs):
    if dryrun:
        print(cmd)
    else:
        try:
            check_call(cmd, shell=True, executable="/usr/bin/zsh", **kwargs)
        except:
            print('##' * 50, '\n', cmd, '\n', "##" * 50)