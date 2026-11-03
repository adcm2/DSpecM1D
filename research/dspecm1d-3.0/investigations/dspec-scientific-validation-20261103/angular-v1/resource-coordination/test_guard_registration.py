"""File-fixture tests only: never start the guardian or control any service."""
from pathlib import Path
import ast
import importlib.util
import json
import tempfile
import unittest

ROOT=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location("private_guard",ROOT/"monitor_tunya_global_ram_v4.py")
guard=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)
def forbidden(*args,**kwargs):
    raise AssertionError("Resource registration tests must not execute process control")
guard.subprocess.run=forbidden

class RegistrationTests(unittest.TestCase):
    def setUp(self):
        self.parent=Path(tempfile.mkdtemp(prefix="dspec-angular-guard-fixture-"))
        self.old=self.parent/"ccdsmextended20261031.slice";self.old.mkdir()
        self.new=self.parent/"dspecscienceangular20261103.slice"
        guard.research_cgroup=lambda:self.old
    def unit(self,root,name,memory=100,freeze="0",pids=2):
        p=root/name;p.mkdir()
        for key,value in {"memory.current":memory,"pids.current":pids,"cgroup.freeze":freeze}.items():
            (p/key).write_text(str(value))
        return p
    def test_old_inventory_unchanged_when_new_slice_absent(self):
        self.unit(self.old,"ccdsm-sp12-owned.service")
        self.unit(self.old,"dspecpilot3d-modal-projection-owned.service")
        self.unit(self.old,"unrelated.service")
        self.unit(self.old,"dspec-science-angular-v1-wrong-root.service")
        self.assertEqual({u["name"] for u in guard.owned_units()},
                         {"ccdsm-sp12-owned.service","dspecpilot3d-modal-projection-owned.service"})
    def test_only_explicit_new_root_and_service_prefix(self):
        self.new.mkdir()
        self.unit(self.new,"dspec-science-angular-v1-run-l800.service",memory=300)
        self.unit(self.new,"ccdsm-sp12-wrong-root.service")
        self.unit(self.new,"unrelated.service")
        self.unit(self.new,"dspec-science-angular-v1-not-service.slice")
        other=self.parent/"unregistered.slice";other.mkdir()
        self.unit(other,"dspec-science-angular-v1-wrong-parent.service")
        units=guard.owned_units()
        self.assertEqual([u["name"] for u in units],["dspec-science-angular-v1-run-l800.service"])
        self.assertEqual(units[0]["path"],str(self.new/units[0]["name"]))
    def test_largest_selection_covers_both_roots(self):
        self.new.mkdir()
        self.unit(self.old,"ccdsm-sp12-small.service",memory=100)
        self.unit(self.new,"dspec-science-angular-v1-large.service",memory=300)
        units=guard.owned_units()
        self.assertEqual(max(units,key=lambda u:u["memory_bytes"])["name"],
                         "dspec-science-angular-v1-large.service")
    def test_missing_original_root_still_fails_closed(self):
        guard.research_cgroup=lambda:self.parent/"missing.slice"
        with self.assertRaises(RuntimeError):guard.owned_units()
    def test_unknown_freezer_state_reported_not_assumed_running(self):
        self.new.mkdir()
        self.unit(self.new,"dspec-science-angular-v1-unknown.service",freeze="invalid")
        units=guard.owned_units()
        self.assertIsNone(units[0]["frozen"])
        self.assertTrue(units[0]["freeze_error"])
    def test_one_pause_and_no_automatic_thaw_policy_unchanged(self):
        self.assertTrue(guard.reconcile_pause(guard.PAUSE_BYTES,None,True,set())["should_pause"])
        self.assertFalse(guard.reconcile_pause(guard.PAUSE_BYTES,None,True,set(),allow_rearm=False)["should_pause"])
        frozen={"dspec-science-angular-v1-running.service"}
        state=guard.reconcile_pause(guard.URGENT_BYTES,next(iter(frozen)),False,frozen)
        self.assertFalse(state["should_pause"])
        self.assertEqual(state["paused_unit"],next(iter(frozen)))
        self.assertTrue(guard.reconcile_pause(guard.PAUSE_BYTES,None,True,{"a","b"})["conflict"])
    def test_only_inventory_function_changed(self):
        manifest=json.loads((ROOT/"supervision-manifest.json").read_text())
        base=Path(manifest["base_source_path"])
        def without_inventory(p):
            tree=ast.parse(p.read_text())
            tree.body=[n for n in tree.body if not(isinstance(n,ast.FunctionDef) and n.name=="owned_units")]
            return ast.dump(tree,include_attributes=False)
        self.assertEqual(without_inventory(base),without_inventory(ROOT/"monitor_tunya_global_ram_v4.py"))

if __name__=="__main__": unittest.main()
